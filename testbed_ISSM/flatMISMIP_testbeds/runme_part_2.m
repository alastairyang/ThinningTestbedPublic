md_idx = [3   4   9  10  15  16]; 
tbl_filename = "runtime_table_2.csv"; 

% parameters
meshsize = 200; % 200 m
perturb_duration = 20; % 20 years
terminus0_x = 56500; % initialized terminus closed to the front
retreat_rate_max = 1000; % maximum retreat rate (m/a)

% read in the table
mdvar_combs = readtable('md_var_combinations.csv');

% initiate a table where we record run time
tableSize = [size(mdvar_combs,1), 4];
varNames = ["geometry", "runtime(min)", "step","finished time"];
varTypes = ["string", "double", "int8","datetime"];
runtimeTbl = table('Size', tableSize, 'VariableTypes', varTypes, 'VariableNames',varNames);

% time
tic
% start iteration
% options: [4,5,6,10,11,12,16,17,18]%[1,2,3,7,8,9,13,14,15]%1:size(mdvar_combs,1)

for jj = md_idx

    var_table = mdvar_combs(jj,:);

    Ly = var_table.('fjord_width');
    delta_gl_depth = var_table.('delta_groundingline_depth');
    bg_fric_coef = var_table.('background_friccoef');

    % names
    modelname = 'yang';
    identifier = ['_W', num2str(Ly), '_GL', num2str(delta_gl_depth), '_FC', num2str(bg_fric_coef)];
    geometry_name = ['domain', identifier];
    par_name = ['par', identifier];

    % make folder for experiment (i.e., this combination of parameters)
    foldername = ['long_models_' modelname '/model' identifier];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end

    % RUN
    for steps = 1:3

        % Cluster parameters
        cluster = generic('name', oshostname(), 'np', 5);
        cluster.interactive = 1;
        waitonlock = 10;

        % Run steps
        org=organizer('repository',foldername,'prefix',['MISMIP_' modelname],'steps',steps);

        if perform(org, 'Mesh_and_Parameterization') % {{{1  STEP 1

            md = model();
            md=parameterize(md, ['long_par_files/' par_name '.par']);
            
            % specify the number of processor based of 1:1000 element ratio
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            cluster.interactive = 1;

            savemodel(org,md);
        end% }}}

        % fixed front initialization
        if perform(org, 'Transient_Steadystate') % {{{1  STEP 2
            % Run the model
            md=loadmodel(org, 'Mesh_and_Parameterization');
        
            md=setflowequation(md,'SSA','all');
            
            % time
            md.timestepping.time_step = 0.010; % default is 0.01
            md.timestepping.final_time=10; % default is 50
            md.settings.output_frequency=100;
        
            md.transient.ismovingfront = 0;
            distance_x = md.mesh.x - terminus0_x;
            if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
            isoline(md,distance_x, 'value', 0, 'output', 'TEMP.exp');
            levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp')); % different from distance_x, levelset is all >0 since it takes abs()
            delete('./TEMP.exp');
            md.mask.ice_levelset(md.mesh.x >terminus0_x) = +abs(levelset(md.mesh.x >terminus0_x));
            md.mask.ice_levelset(md.mesh.x<=terminus0_x) = -abs(levelset(md.mesh.x<=terminus0_x));
        
            md.stressbalance.abstol=NaN;
            md.stressbalance.restol=1e-4;
            md.stressbalance.maxiter=50;
            md.verbose=verbose('convergence',false,'solution',true);
            md.cluster=cluster;
            md.toolkits.DefaultAnalysis = bcgslbjacobioptions(); % faster solution; robust for 2D SSA
            md.settings.solver_residue_threshold = 1e-5; % I set it; Denis didn't.
            md.settings.waitonlock=waitonlock;

            % config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;

            % solve
            md=solve(md,'tr');
            if md.settings.waitonlock == 0
                fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
                fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
                return
            end
            
            % extend, but change the boundary U and H combinations
            % to mitigate the shock near the influx boundary
            md = transientrestart(md);

            % parameters for the grid
            ds = 50;
            Lx = max(md.mesh.x);
            Ly = max(md.mesh.y);
            x = 0:ds:Lx-ds;
            y = 0:ds:Ly-ds;
            H = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                     md.geometry.thickness,...
                                     x, y, NaN);
            U = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                     md.initialization.vx,...
                                     x, y, NaN);
            % get the velocity and thickness
            dist = 2000; % a semi-arbitrary distance downstream of where shock usually is
            h_slice = H(:,dist/ds);
            u_slice = U(:,dist/ds);
            % find the boundary elements
            pos=find(md.mesh.x<0.1 & md.mesh.x>-0.1);
            y_axis = md.mesh.y(pos);
            % interpolate flux, h and get the corresponding v
            uh_x0 = interp1(y,u_slice.*h_slice, y_axis,'spline');
            h_x0 = interp1(y,h_slice,y_axis,'spline');
            v_x0 = uh_x0./h_x0;
            
            % assign to md
            md.masstransport.spcthickness(pos) = h_x0;
            md.stressbalance.spcvx(pos) = v_x0;
            md.stressbalance.spcvx(md.stressbalance.spcvx<0) = 0;
            
            % cluster
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            
            % time
            md.settings.output_frequency = 50;
            md.timestepping.final_time = md.timestepping.start_time + 10;
            
            % solve
            md = solve(md,'tr');

            savemodel(org,md);
    
            % write time
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end% }}}

        if perform(org, 'Transient_Steadystate_Extended') % {{{1  STEP 3
            % Extend the initialization runs with adaptive time stepping.
            md = loadmodel(org, 'Transient_Steadystate');
            md = transientrestart(md);

            % switch to adaptive time stepping and 
            start_time = md.timestepping.start_time;
            md.timestepping = timesteppingadaptive();
            md.timestepping.time_step_min = 0.01;
            md.timestepping.start_time = start_time;
            md.timestepping.final_time = md.timestepping.start_time + 100;
            
            % parameters
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            
            % time
            md.settings.output_frequency = 100;
            md.transient.requested_outputs={'default','IceVolume'};

            % solve and save
            md = solve(md,'tr');
            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end% }}}

        if perform(org, 'Transient_ExtraInfo') % {{{1  STEP 4
        % this step merely acquires some stress and strain rate data for
        % the next steps.
            md = loadmodel(org, 'Transient_Steadystate_Extended');
            md = transientrestart(md);

            % extend the run time for one year
            md.timestepping.final_time = md.timestepping.start_time + 1;
            md.settings.output_frequency = 1;

            % output
            md.stressbalance.requested_outputs={'default','StrainRatexx',...
                                                'StrainRateyy','StrainRatexy',...
                                                'StrainRateeffective', ...
                                                'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            md = solve(md,'tr');

            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end

        if perform(org, 'Transient_CalvingOnly')% {{{1 STEP 5
            md = loadmodel(org, 'Transient_ExtraInfo');
            md = transientrestart(md);

            % revert back to constant time stepping (for easier data
            % analysis later / ensure stability)
            end_time = perturb_duration;
            start_time = md.timestepping.start_time;
            md.timestepping = timestepping(); 
            md.timestepping.time_step = 0.01;
            md.timestepping.start_time = start_time + md.timestepping.time_step; % displace by 1 dt, to align with calving_shearmargin sim results
            md.timestepping.final_time = start_time + end_time;
            md.settings.output_frequency = 10;
            
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_sequence = [retreat_advance, retreat_slow];

            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), even if we are
            % prescribing the terminus
            md.transient.ismovingfront = 1;
            md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);

            % create sequences of terminus position via spclevelset
            levelset0 = md.mask.ice_levelset;
            md.levelset.spclevelset = [];
            md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
            
            culmu_magnitude = 0;
            calving_start = 1;
            calving_end   = end_time;
            
            % the year before it starts: still at the original position 
            %md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time + calving_start - 1];
            for time = calving_start : calving_end % only the few years in the middle
                magnitude = retreat_sequence(time);
                culmu_magnitude = culmu_magnitude + magnitude;
                signeddistance = move_terminus_levelset_mod(md, levelset0, culmu_magnitude, -1, true);

                signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
                pos      = find(signeddistance<0);

                if exist('TEMP.exp','file'), delete('TEMP.exp'); end
                isoline(md, signeddistance, 'value', 0, 'output', 'TEMP.exp');
                signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
                delete('TEMP.exp');
                signeddistance(pos) = -signeddistance(pos);

                md.levelset.spclevelset(:,end+1) = [signeddistance; md.timestepping.start_time + time];
            end
            
            % cluster config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % solve
            md = solve(md,'tr');

            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])

        end%}}}
        
        if perform(org, 'Transient_Calving_SMweakening')% {{{1 STEP 6
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_smw = 0.1; % shear margin weakening update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_sequence = [retreat_advance, retreat_slow];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), even if we are
            % prescribing the terminus
            md.transient.ismovingfront = 1;
            md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);
            
            % create sequences of terminus position via spclevelset
            levelset0 = md.mask.ice_levelset;
            md.levelset.spclevelset = [];
            md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
            
            culmu_magnitude = 0;
            calving_start = 1;
            calving_end   = end_time;
            
            % prescribing calving front: change annually
            for time = calving_start : dt_calve: calving_end % only the few years in the middle
                magnitude = retreat_sequence(time);
                culmu_magnitude = culmu_magnitude + magnitude;
                signeddistance = move_terminus_levelset_mod(md, levelset0, culmu_magnitude, -1, true);

                signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
                pos = find(signeddistance<0);

                if exist('TEMP.exp','file'), delete('TEMP.exp'); end
                isoline(md, signeddistance, 'value', 0, 'output', 'TEMP.exp');
                signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
                delete('TEMP.exp');
                signeddistance(pos) = -signeddistance(pos);

                md.levelset.spclevelset(:,end+1) = [signeddistance; md.timestepping.start_time + time];
            end
            
            %% Shear margin weakening

            % initialize B, T, effective strain rate (ee)
            T = 273.15*ones(md.mesh.numberofvertices, 1); % initialize T (melting point)
            B0 = md.materials.rheology_B;
            B = [md.materials.rheology_B;md.timestepping.start_time];
            md.materials.rheology_B = B;
            ee0 = md.results.TransientSolution(end).StrainRateeffective;

            % save previous fields separately
            % this step help re-assembles all results later easily
            md_temp = transientrestart(md);
            previous_results = md_temp.results;
            clear md_temp
            
            % initialize
            new_results = [];
            dB = zeros(size(ee0));

            for it = 1:end_time/dt_smw                                    
                ee = md.results.TransientSolution(end).StrainRateeffective;

                % restart and specify sim duration
                md = transientrestart(md);
                md.timestepping.time_step = 0.01;
                md.timestepping.final_time = md.timestepping.start_time + dt_smw;
                md.settings.output_frequency = dt_smw/md.timestepping.time_step;

                if it == 1 % first iteration: do not update B. Just run for another dt
                    current_time = md.timestepping.start_time;
                    B = B(1:end-1);
                    B = [B; current_time + dt_smw];
                    md.materials.rheology_B = [md.materials.rheology_B, B];

                    md = solve(md,'tr');
                    % save
                    new_results = [new_results,md.results.TransientSolution(1:end-1)];
                else
                    % calculate new B and update the model
                    %dee = ee*factor - ee0*factor; 
                    B = B(1:end-1);
                    ice_mask = md.mask.ice_levelset;
                    [B, T, ~, dB] = softening_map_studies( B, ee, T, dt_smw, ee0, B0, ice_mask);
                    % append time and assign
                    current_time = md.timestepping.start_time;
                    B = [B; current_time + dt_smw];
                    md.materials.rheology_B = [md.materials.rheology_B, B];
                    
                    % solve
                    md = solve(md,'tr');
    
                    % save the new result to a separate var
                    new_results = [new_results,md.results.TransientSolution(1)];
                end
            end
            md.results = previous_results;
            md.results.TransientSolution = new_results;
            clear previous_results

            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end

        if perform(org, 'Transient_Calving_MassUnloading')% {{{1 STEP 7
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_sequence = [retreat_advance, retreat_slow];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), even if we are
            % prescribing the terminus
            md.transient.ismovingfront = 1;
            md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);
            
            % create sequences of terminus position via spclevelset
            levelset0 = md.mask.ice_levelset;
            md.levelset.spclevelset = [];
            md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
            
            culmu_magnitude = 0;
            calving_start = 1;
            calving_end   = end_time;
            
            % prescribing calving front: change annually
            for time = calving_start : dt_calve: calving_end % only the few years in the middle
                magnitude = retreat_sequence(time);
                culmu_magnitude = culmu_magnitude + magnitude;
                signeddistance = move_terminus_levelset_mod(md, levelset0, culmu_magnitude, -1, true);

                signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
                pos = find(signeddistance<0);

                if exist('TEMP.exp','file'), delete('TEMP.exp'); end
                isoline(md, signeddistance, 'value', 0, 'output', 'TEMP.exp');
                signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
                delete('TEMP.exp');
                signeddistance(pos) = -signeddistance(pos);

                md.levelset.spclevelset(:,end+1) = [signeddistance; md.timestepping.start_time + time];
            end
            
            %% Mass unloading
            % save previous fields separately
            % this step help re-assembles all results later easily
            md_temp = transientrestart(md);
            previous_results = md_temp.results;
            next_start_time = md_temp.timestepping.start_time;
            clear md_temp
            
            % initialize
            new_results = [];

            % get the equivalent coefficients if using Budd sliding law
            law_from = 'Weertman';
            law_to = 'Budd';
            C0 = md.friction.C;
            H0 = md.results.TransientSolution(end).Thickness;
            Zb = md.results.TransientSolution(end).Base;
            k_budd = fric_coef_conversion(law_from, law_to, md, C0, H0, Zb);

            % add an initial time to the friction coef vector
            md.friction.C = [C0; next_start_time];

            for it = 1:end_time/dt_mu                                   
                results = md.results.TransientSolution;
                % restart and specify sim duration
                md = transientrestart(md);
                md.timestepping.time_step = 0.01;
                md.timestepping.final_time = md.timestepping.start_time + dt_mu;
                md.settings.output_frequency = dt_mu/md.timestepping.time_step;

                % calculate new fric coef
                if it == 1
                    deltaH = zeros(size(md.geometry.thickness));
                    C = C0;
                else
                    deltaH = results(end).Thickness - H0;
                end
                ocean_mask = results(end).MaskOceanLevelset;
                C = mass_unloading(md, deltaH, k_budd, C0, C, ocean_mask);
                % append time and assign
                current_time = md.timestepping.start_time;
                C_add_time = [C; current_time + dt_mu];
                md.friction.C = [md.friction.C, C_add_time];

                % solve
                md = solve(md,'tr');

                % save the new result to a separate var
                new_results = [new_results,md.results.TransientSolution(1)];
            end
            md.results = previous_results;
            md.results.TransientSolution = new_results;
            clear previous_results new_results

            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            runtimeTbl{jj,1} = string(geometry_name);
            runtimeTbl{jj,2} = runTime/60;
            runtimeTbl{jj,3} = steps;
            runtimeTbl{jj,4} = datetime;
            writetable(runtimeTbl, tbl_filename);
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end
           

    end
end
