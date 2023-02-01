% runme_free_front.m template.
%%%%%%%%%%%%%% DON'T CHANGE THE FOLLOWING LINES / DON'T MODIFY LINE NO. %%%%%%%%%%%%%%
% specify model index; specify 
md_idx = 1:18;
tbl_filename = "runtime_table_all.csv";
gauss_mag = 0.8; % fractional reduction in basal shear stress, (1-gauss_mag)*tau_b
%%%%%%%%%%%%%% 

% parameters
meshsize = 200; % size of mesh element (m)
perturb_duration = 16; % 16 years; but it's not the total perturbation run duration
no_retreat_duration = 5; % 5 years no retreat, in the beginning & end
terminus0_x = 56500; % initialized terminus closed to the front
retreat_rate_max = 1000; % maximum retreat rate (m/a)
seasonal_retreat_rate_max = 120; % maximum retreat rate for seasonal calving experiment (m/a)
calve_seasonal_max = 50; % seasonal variation 
gauss_xloc = 3.2e4; % x-axis location of the gaussian perturbation (m)
%gauss_tscale = 2; % flucation time scale. One full cyle is 2*2 years.
gauss_width_ratio = 0.08; % the ratio of gaussian perturbation patch to fjord width
gauss_timestep = 0.01; % temporal resolution of the prescribed gaussian path (yr)
pulse_gauss_tscale = 0.1; % timescale of the pulse gaussian perturbation (yr)
diffu_gauss_tscale = 1.3; % timescale of the diffused gaussian perturbation (yr)
gauss_efold = 3; % e-folding, decay to e^(-n) 
gauss_perturb_repeat_tscale = 2; % repeat time period for both pulse and diffused gaussian perturb
pulse_gauss_tshift = 0.5; % for each repeat t cycle, the pulse is placed at the given time (yr) after the start
ds = 100; % grid spacing (m), when interpolate from mesh to grid

% Gaussian basal perturbation simulation titles
pulse_gauss_title = ['Transient_Calving_PulseGaussianPerturb_', num2str(gauss_mag*10)];
diffu_gauss_title = ['Transient_Calving_DiffuGaussianPerturb_', num2str(gauss_mag*10)];
pulse_gauss_mu_title = ['Transient_Calving_MassUnloading_PulseGaussianPerturb_', num2str(gauss_mag*10)];
diffu_gauss_mu_title = ['Transient_Calving_MassUnloading_DiffuGaussianPerturb_', num2str(gauss_mag*10)];

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
    for steps = 11

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
            md.timestepping.final_time = md.timestepping.start_time + 500;
            
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
            end_time = perturb_duration + 2*no_retreat_duration;
            start_time = md.timestepping.start_time;
            md.timestepping = timestepping(); 
            md.timestepping.time_step = 0.01;
            md.timestepping.start_time = start_time + md.timestepping.time_step; % displace by 1 dt, to align with calving_shearmargin sim results
            md.timestepping.final_time = start_time + end_time;
            md.settings.output_frequency = 10;
            
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];

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
            for time = calving_start : calving_end
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
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_smw = 0.1; % shear margin weakening update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
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
                    try
                        md = solve(md,'tr');
                    catch
                        pause
                    end
    
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
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
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
            % mass unloading activation time
            % we only allow this effective pressure feedback (mass
            % unloading) to be active after the terminus retreat has
            % started. 
            mu_time_mask = zeros(size(retreat_sequence));
            mu_time_mask(find(retreat_sequence > 0, 1,'first'):end) = 1;
            % interp
            mu_time_mask_interp = interp1(1:end_time, mu_time_mask, 0:dt_mu:end_time-dt_mu, 'previous',0); 
        
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
                if it == 1 % initial condition: delta(H) = 0
                    deltaH = mu_time_mask_interp(it)*zeros(size(md.geometry.thickness));
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0);
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

        if perform(org, 'Transient_SeasonalCalving_MassUnloading')% {{{1 STEP 8
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 0.1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(0,seasonal_retreat_rate_max, perturb_duration/2*(1/dt_calve));
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration*(1/dt_calve));
            retreat_perturb = [retreat_advance, retreat_slow];
            retreat_seasonal = calve_seasonal_max*sqrt(retreat_perturb).*sin(2*pi*[0:0.1:perturb_duration-0.1]);
            retreat_sequence = [retreat_no, retreat_perturb + retreat_seasonal, retreat_no];
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
            for time_i = 1:length(retreat_sequence) % only the few years in the middle
                magnitude = retreat_sequence(time_i);
                culmu_magnitude = culmu_magnitude + magnitude;
                signeddistance = move_terminus_levelset_mod(md, levelset0, culmu_magnitude, -1, true);

                signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
                pos = find(signeddistance<0);

                if exist('TEMP.exp','file'), delete('TEMP.exp'); end
                isoline(md, signeddistance, 'value', 0, 'output', 'TEMP.exp');
                signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
                delete('TEMP.exp');
                signeddistance(pos) = -signeddistance(pos);

                md.levelset.spclevelset(:,end+1) = [signeddistance; md.timestepping.start_time + time_i*dt_calve];
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

            % mass unloading activation time
            % we only allow this effective pressure feedback (mass
            % unloading) to be active after the terminus retreat has
            % started. 
            mu_time_mask = zeros(size(retreat_sequence));
            mu_time_mask(find(retreat_sequence > 0, 1,'first'):end) = 1;
            % interp
            mu_time_mask_interp = interp1(1:end_time, mu_time_mask, 0:dt_mu:end_time-dt_mu, 'previous',0); 

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
                if it == 1 % initial condition: delta(H) = 0
                    deltaH = mu_time_mask_interp(it)*zeros(size(md.geometry.thickness));
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0);
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

        if perform(org, pulse_gauss_title)% {{{1 STEP 9: Basal perturbation: pulse-like gaussian patch
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), as we are
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
            
            %% Basal perturbation with a Gaussian patch
            % add an initial time column to the friction coef vector
            C0 = md.friction.C;
            md.friction.C = [C0; start_time];
            init_taub = C0.^2.*md.results.TransientSolution(end).Vel./md.constants.yts;
            % create the temporal fluctuation sequence, delta_tau
            start_t = md.levelset.spclevelset(end,1);
            perturb_t = 0:gauss_timestep:perturb_duration-gauss_timestep;
            no_perturb_t = 0:gauss_timestep:no_retreat_duration-gauss_timestep;

            % get pulse gaussian (time dimension)
            pulse_gauss = gauss_mag*make_pulse_gauss(pulse_gauss_tscale, gauss_efold, gauss_perturb_repeat_tscale, gauss_timestep, pulse_gauss_tshift);
            % stack them to make a full sequence
            total_cycle = perturb_duration/gauss_perturb_repeat_tscale;
            gauss_mags = repmat(pulse_gauss, 1, total_cycle);
            % add the no perturb period
            gauss_mags = [zeros(size(no_perturb_t)), gauss_mags, zeros(size(no_perturb_t))];
            % actual time axis
            gauss_t = 0:gauss_timestep:(perturb_duration+2*no_retreat_duration-gauss_timestep);
            gauss_t = gauss_t + (start_t + gauss_timestep);
            % interp from mesh to grid and smooth
            Lx = max(md.mesh.x);
            Ly = max(md.mesh.y);
            x = 0:ds:Lx;
            y = 0:ds:Ly;
            [X,Y] = meshgrid(x, y);
            % find centerline index
            if rem(size(X,1), 2) == 0
                mid_i = size(X,1)/2;
            else
                mid_i = (size(X,1)+1)/2;
            end
            init_taub_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            init_taub,x, y, NaN);
            % smooth the initial basal shear stress field
            for row = 1:size(init_taub_grid,1)
                init_taub_grid(row,:) = smooth(init_taub_grid(row,:),20);
            end
            % location of the the perturbation
            x0 = gauss_xloc;
            y0_i = mid_i;
            x0_i = x0/ds;
            y0 = ds*mid_i;
            % Friction coefficient fields for all time steps
            % record the results at every 0.1 yr
            for iter = 1:length(gauss_t)
                amp = gauss_mags(iter)*init_taub_grid(y0_i,x0_i);
                % scale with respect to the max width in our testbeds
                max_width = max(mdvar_combs.fjord_width);
                width_ratio = var_table.('fjord_width')/max_width;
                width = gauss_width_ratio*max_width*sqrt(width_ratio);
                delta_taub = transient_slippatch(X,Y,x0,y0,width,amp);
                % convert back to changes in fric coefficient
                delta_taub = InterpFromGridToMesh(x',y',delta_taub,md.mesh.x,md.mesh.y,0);
                delta_C = sqrt(delta_taub./(md.results.TransientSolution(end).Vel/md.constants.yts));
                delta_C(isinf(delta_C)) = 0;
                delta_C(isnan(delta_C)) = 0;
                C_new = [C0 - delta_C; gauss_t(iter)];
                md.friction.C = [md.friction.C, C_new];
            end

            % restart and specify sim duration
            md = transientrestart(md);
            md.timestepping.time_step = 0.01;
            md.timestepping.final_time = md.timestepping.start_time + end_time;
            md.settings.output_frequency = 10;

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
        end

        if perform(org, diffu_gauss_title)% {{{1 STEP 10: Basal perturbation: diffused gaussian patch
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), as we are
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
            
            %% Basal perturbation with a Gaussian patch
            % add an initial time column to the friction coef vector
            C0 = md.friction.C;
            md.friction.C = [C0; start_time];
            init_taub = C0.^2.*md.results.TransientSolution(end).Vel./md.constants.yts;
            % create the temporal fluctuation sequence, delta_tau
            start_t = md.levelset.spclevelset(end,1);
            perturb_t = 0:gauss_timestep:perturb_duration-gauss_timestep;
            no_perturb_t = 0:gauss_timestep:no_retreat_duration-gauss_timestep;

            % get diffused pulse gaussian
            diffu_gauss = gauss_mag*make_diffu_gauss(pulse_gauss_tscale, diffu_gauss_tscale, gauss_efold, gauss_perturb_repeat_tscale, gauss_timestep);
            % stack them to make a full sequence
            total_cycle = perturb_duration/gauss_perturb_repeat_tscale;
            gauss_mags = repmat(diffu_gauss, 1, total_cycle);
            % add the no perturb period
            gauss_mags = [zeros(size(no_perturb_t)), gauss_mags, zeros(size(no_perturb_t))];
            % actual time axis
            gauss_t = 0:gauss_timestep:(perturb_duration+2*no_retreat_duration-gauss_timestep);
            gauss_t = gauss_t + (start_t + gauss_timestep);
            % interp from mesh to grid and smooth
            Lx = max(md.mesh.x);
            Ly = max(md.mesh.y);
            x = 0:ds:Lx;
            y = 0:ds:Ly;
            [X,Y] = meshgrid(x, y);
            % find centerline index
            if rem(size(X,1), 2) == 0
                mid_i = size(X,1)/2;
            else
                mid_i = (size(X,1)+1)/2;
            end
            init_taub_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            init_taub,x, y, NaN);
            % smooth the initial basal shear stress field
            for row = 1:size(init_taub_grid,1)
                init_taub_grid(row,:) = smooth(init_taub_grid(row,:),20);
            end
            % location of the the perturbation
            x0 = gauss_xloc;
            y0_i = mid_i;
            x0_i = x0/ds;
            y0 = ds*mid_i;
            % Friction coefficient fields for all time steps
            % record the results at every 0.1 yr
            for iter = 1:length(gauss_t)
                amp = gauss_mags(iter)*init_taub_grid(y0_i,x0_i);
                % scale with respect to the max width in our testbeds
                max_width = max(mdvar_combs.fjord_width);
                width_ratio = var_table.('fjord_width')/max_width;
                width = gauss_width_ratio*max_width*sqrt(width_ratio);
                delta_taub = transient_slippatch(X,Y,x0,y0,width,amp);
                % convert back to changes in fric coefficient
                delta_taub = InterpFromGridToMesh(x',y',delta_taub,md.mesh.x,md.mesh.y,0);
                delta_C = sqrt(delta_taub./(md.results.TransientSolution(end).Vel/md.constants.yts));
                delta_C(isinf(delta_C)) = 0;
                delta_C(isnan(delta_C)) = 0;
                C_new = [C0 - delta_C; gauss_t(iter)];
                md.friction.C = [md.friction.C, C_new];
            end

            % restart and specify sim duration
            md = transientrestart(md);
            md.timestepping.time_step = 0.01;
            md.timestepping.final_time = md.timestepping.start_time + end_time;
            md.settings.output_frequency = 10;

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
        end

        if perform(org, pulse_gauss_mu_title)% {{{1 STEP 11: Basal perturbation: pulse gaussian patch + mass unloading
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), as we are
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
            
            %% Basal perturbation with a Gaussian patch
            % add an initial time column to the friction coef vector
            C0 = md.friction.C;
            init_taub = C0.^2.*md.results.TransientSolution(end).Vel./md.constants.yts;
            % create the temporal fluctuation sequence, delta_tau
            start_t = md.levelset.spclevelset(end,1);
            perturb_t = 0:gauss_timestep:perturb_duration-gauss_timestep;
            no_perturb_t = 0:gauss_timestep:no_retreat_duration-gauss_timestep;

            % get pulse gaussian
            pulse_gauss = gauss_mag*make_pulse_gauss(pulse_gauss_tscale, gauss_efold, gauss_perturb_repeat_tscale, gauss_timestep, pulse_gauss_tshift);
            % stack them to make a full sequence
            total_cycle = perturb_duration/gauss_perturb_repeat_tscale;
            gauss_mags = repmat(pulse_gauss, 1, total_cycle);
            % add the no perturb period
            gauss_mags = [zeros(size(no_perturb_t)), gauss_mags, zeros(size(no_perturb_t))];
            % actual time axis
            gauss_t = 0:gauss_timestep:(perturb_duration+2*no_retreat_duration-gauss_timestep);
            gauss_t = gauss_t + (start_t + gauss_timestep);
            % interp from mesh to grid and smooth
            Lx = max(md.mesh.x);
            Ly = max(md.mesh.y);
            x = 0:ds:Lx;
            y = 0:ds:Ly;
            [X,Y] = meshgrid(x, y);
            % find centerline index
            if rem(size(X,1), 2) == 0
                mid_i = size(X,1)/2;
            else
                mid_i = (size(X,1)+1)/2;
            end
            init_taub_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            init_taub,x, y, NaN);
            % smooth the initial basal shear stress field
            for row = 1:size(init_taub_grid,1)
                init_taub_grid(row,:) = smooth(init_taub_grid(row,:),20);
            end
            % location of the the perturbation
            x0 = gauss_xloc;
            y0_i = mid_i;
            x0_i = x0/ds;
            y0 = ds*mid_i;
            % Friction coefficient fields for all time steps
            % save to a matrix
            delta_C_all = zeros(length(C0), length(gauss_t));
            for iter = 1:length(gauss_t)
                amp = gauss_mags(iter)*init_taub_grid(y0_i,x0_i);
                % scale with respect to the max width in our testbeds
                max_width = max(mdvar_combs.fjord_width);
                width_ratio = var_table.('fjord_width')/max_width;
                width = gauss_width_ratio*max_width*sqrt(width_ratio);
                delta_taub = transient_slippatch(X,Y,x0,y0,width,amp);
                % convert back to changes in fric coefficient
                delta_taub = InterpFromGridToMesh(x',y',delta_taub,md.mesh.x,md.mesh.y,0);
                delta_C = sqrt(delta_taub./(md.results.TransientSolution(end).Vel/md.constants.yts));
                delta_C(isinf(delta_C)) = 0;
                delta_C(isnan(delta_C)) = 0;
                delta_C_all(:,iter) = delta_C; % delta_C: reduction is positive
                % this delta_C_all is then modified below 
                % when we calculate the effective pressure change
                % and then assimilate into md.friction.C
            end

            %% Mass unloading
            % save previous fields separately
            % this step help re-assembles all results later easily
            md_temp = transientrestart(md);
            previous_results = md_temp.results;
            next_start_time = md_temp.timestepping.start_time;
            clear md_temp
            
            % initialize a var to collect results
            new_results = [];
            % mass unloading activation time
            % we only allow this effective pressure feedback (mass
            % unloading) to be active after the terminus retreat has
            % started. 
            mu_time_mask = zeros(size(retreat_sequence));
            mu_time_mask(find(retreat_sequence > 0, 1,'first'):end) = 1;
            % interp
            mu_time_mask_interp = interp1(1:end_time, mu_time_mask, 0:gauss_timestep:end_time-gauss_timestep, 'previous',0); 
        
            % get the equivalent coefficients if using Budd sliding law
            law_from = 'Weertman';
            law_to = 'Budd';
            C0 = md.friction.C;
            H0 = md.results.TransientSolution(end).Thickness;
            Zb = md.results.TransientSolution(end).Base;
            k_budd = fric_coef_conversion(law_from, law_to, md, C0, H0, Zb);

            % add an initial time to the friction coef vector
            md.friction.C = [C0; next_start_time];

            it_count = 0;
            for it = 1:(end_time/gauss_timestep)-1
                it_count = it_count + 1;
                results = md.results.TransientSolution;
                % restart and specify sim duration
                md = transientrestart(md);
                md.timestepping.time_step = 0.01;
                md.timestepping.final_time = md.timestepping.start_time + gauss_timestep;
                md.settings.output_frequency = 1;

                % calculate new fric coef
                if it == 1 % initial condition: delta(H) = 0
                    deltaH = mu_time_mask_interp(it)*zeros(size(md.geometry.thickness));
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0);
                end
                ocean_mask = results(end).MaskOceanLevelset;
                C = mass_unloading(md, deltaH, k_budd, C0, C, ocean_mask);
                C = C - delta_C_all(:,it);
                C(C<0) = 0;
                % append time and assign
                current_time = md.timestepping.start_time;
                C_add_time = [C; current_time + dt_mu];
                md.friction.C = [md.friction.C, C_add_time];

                % solve
                md = solve(md,'tr');
                
                if it_count == 10 || it == 1 % save every 10 time steps (every 0.1 year) & save the first time
                    it_count = 0;
                    new_results = [new_results, md.results.TransientSolution(1)];
                    % remove the extra fields
                    names = fieldnames(md.results);
                    results_cell = struct2cell(md.results);
                    md.results = cell2struct(results_cell(names == "TransientSolution"),...
                                             names(names == "TransientSolution"));
                    clear names results_cell
                end
            end
            md.results = previous_results;
            md.results.TransientSolution = new_results;

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

        if perform(org, diffu_gauss_mu_title)% {{{1 STEP 12: Basal perturbation: diffused pulse gaussian patch + mass unloading
            md = loadmodel(org, 'Transient_ExtraInfo');

            % parameter regarding time
            end_time = perturb_duration + 2*no_retreat_duration;

            start_time = md.timestepping.final_time;
            md.timestepping = timestepping(); 
            md.timestepping.start_time = start_time;
            dt_mu = 0.1; % mass unloading update dt
            dt_calve = 1; % calving front position update dt

            % simulation config
            np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
            cluster = generic('name', oshostname(), 'np', np);
            md.cluster = cluster;
            % relax max iteration (might need in certain shear margin runs)
            md.stressbalance.maxiter=100;
            % do not interpolate forcing
            md.timestepping.interp_forcing = 0;

            %% Calving
            % forcings
            retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
            retreat_slow = flip(retreat_advance);
            retreat_no = zeros(1,no_retreat_duration);
            retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);

            % enabling movingfront (levelset method), as we are
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
            
            %% Basal perturbation with a Gaussian patch
            % add an initial time column to the friction coef vector
            C0 = md.friction.C;
            init_taub = C0.^2.*md.results.TransientSolution(end).Vel./md.constants.yts;
            % create the temporal fluctuation sequence, delta_tau
            start_t = md.levelset.spclevelset(end,1);
            perturb_t = 0:gauss_timestep:perturb_duration-gauss_timestep;
            no_perturb_t = 0:gauss_timestep:no_retreat_duration-gauss_timestep;

            % get diffused pulse gaussian
            diffu_gauss = gauss_mag*make_diffu_gauss(pulse_gauss_tscale, diffu_gauss_tscale, gauss_efold, gauss_perturb_repeat_tscale, gauss_timestep);
            % stack them to make a full sequence
            total_cycle = perturb_duration/gauss_perturb_repeat_tscale;
            gauss_mags = repmat(diffu_gauss, 1, total_cycle);
            % add the no perturb period
            gauss_mags = [zeros(size(no_perturb_t)), gauss_mags, zeros(size(no_perturb_t))];
            % actual time axis
            gauss_t = 0:gauss_timestep:(perturb_duration+2*no_retreat_duration-gauss_timestep);
            gauss_t = gauss_t + (start_t + gauss_timestep);
            % interp from mesh to grid and smooth
            Lx = max(md.mesh.x);
            Ly = max(md.mesh.y);
            x = 0:ds:Lx;
            y = 0:ds:Ly;
            [X,Y] = meshgrid(x, y);
            % find centerline index
            if rem(size(X,1), 2) == 0
                mid_i = size(X,1)/2;
            else
                mid_i = (size(X,1)+1)/2;
            end
            init_taub_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            init_taub,x, y, NaN);
            % smooth the initial basal shear stress field
            for row = 1:size(init_taub_grid,1)
                init_taub_grid(row,:) = smooth(init_taub_grid(row,:),20);
            end
            % location of the the perturbation
            x0 = gauss_xloc;
            y0_i = mid_i;
            x0_i = x0/ds;
            y0 = ds*mid_i;
            % Friction coefficient fields for all time steps
            % save to a matrix
            delta_C_all = zeros(length(C0), length(gauss_t));
            for iter = 1:length(gauss_t)
                amp = gauss_mags(iter)*init_taub_grid(y0_i,x0_i);
                % scale with respect to the max width in our testbeds
                max_width = max(mdvar_combs.fjord_width);
                width_ratio = var_table.('fjord_width')/max_width;
                width = gauss_width_ratio*max_width*sqrt(width_ratio);
                delta_taub = transient_slippatch(X,Y,x0,y0,width,amp);
                % convert back to changes in fric coefficient
                delta_taub = InterpFromGridToMesh(x',y',delta_taub,md.mesh.x,md.mesh.y,0);
                delta_C = sqrt(delta_taub./(md.results.TransientSolution(end).Vel/md.constants.yts));
                delta_C(isinf(delta_C)) = 0;
                delta_C(isnan(delta_C)) = 0;
                delta_C_all(:,iter) = delta_C; % delta_C: reduction is positive
                % this delta_C_all is then modified below 
                % when we calculate the effective pressure change
                % and then assimilate into md.friction.C
            end

            %% Mass unloading
            % save previous fields separately
            % this step help re-assembles all results later easily
            md_temp = transientrestart(md);
            previous_results = md_temp.results;
            next_start_time = md_temp.timestepping.start_time;
            clear md_temp
            
            % initialize a var to collect results
            new_results = [];
            % mass unloading activation time
            % we only allow this effective pressure feedback (mass
            % unloading) to be active after the terminus retreat has
            % started. 
            mu_time_mask = zeros(size(retreat_sequence));
            mu_time_mask(find(retreat_sequence > 0, 1,'first'):end) = 1;
            % interp
            mu_time_mask_interp = interp1(1:end_time, mu_time_mask, 0:gauss_timestep:end_time-gauss_timestep, 'previous',0); 
        
            % get the equivalent coefficients if using Budd sliding law
            law_from = 'Weertman';
            law_to = 'Budd';
            C0 = md.friction.C;
            H0 = md.results.TransientSolution(end).Thickness;
            Zb = md.results.TransientSolution(end).Base;
            k_budd = fric_coef_conversion(law_from, law_to, md, C0, H0, Zb);

            % add an initial time to the friction coef vector
            md.friction.C = [C0; next_start_time];

            it_count = 0;
            for it = 1:(end_time/gauss_timestep)-1
                it_count = it_count + 1;
                results = md.results.TransientSolution;
                % restart and specify sim duration
                md = transientrestart(md);
                md.timestepping.time_step = 0.01;
                md.timestepping.final_time = md.timestepping.start_time + gauss_timestep;
                md.settings.output_frequency = 1;

                % calculate new fric coef
                if it == 1 % initial condition: delta(H) = 0
                    deltaH = mu_time_mask_interp(it)*zeros(size(md.geometry.thickness));
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0);
                end
                ocean_mask = results(end).MaskOceanLevelset;
                C = mass_unloading(md, deltaH, k_budd, C0, C, ocean_mask);
                C = C - delta_C_all(:,it);
                C(C<0) = 0;
                % append time and assign
                current_time = md.timestepping.start_time;
                C_add_time = [C; current_time + dt_mu];
                md.friction.C = [md.friction.C, C_add_time];

                % solve
                md = solve(md,'tr');
                
                if it_count == 10 || it == 1 % save every 10 time steps (every 0.1 year)
                    it_count = 0;
                    new_results = [new_results, md.results.TransientSolution(1)];
                    % remove the extra fields
                    names = fieldnames(md.results);
                    results_cell = struct2cell(md.results);
                    md.results = cell2struct(results_cell(names == "TransientSolution"),...
                                             names(names == "TransientSolution"));
                    clear names results_cell
                end
            end
            md.results = previous_results;
            md.results.TransientSolution = new_results;

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



