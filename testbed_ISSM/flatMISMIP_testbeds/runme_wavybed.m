%%%%%%%%%%%%%% DON'T CHANGE THE FOLLOWING LINES / DON'T MODIFY LINE NO. %%%%%%%%%%%%%%
% specify model index; specify 
md_idx = [6,18];
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

table_vals  = [meshsize,  perturb_duration, no_retreat_duration,   terminus0_x,  retreat_rate_max,...
              seasonal_retreat_rate_max, calve_seasonal_max,   gauss_xloc, gauss_mag, gauss_width_ratio, ...
              gauss_timestep, pulse_gauss_tscale,   diffu_gauss_tscale,  gauss_efold, gauss_perturb_repeat_tscale,...
              pulse_gauss_tshift, ds];
table_names = ["meshsize","perturb_duration","no_retreat_duration","terminus0_x","retreat_rate_max",...
             "seasonal_retreat_rate_max","calve_seasonal_max","gauss_xloc","gauss_mag","gauss_width_ratio",...
             "gauss_timestep","pulse_gauss_tscale","diffu_gauss_tscale","gauss_efold", "gauss_perturb_repeat_tscale",...
             "pulse_gauss_tshift","ds"];
runme_param = array2table(table_vals, 'VariableNames',table_names);
writetable(runme_param, 'runme_param.csv')

% Gaussian basal perturbation simulation titles
pulse_gauss_title = ['Transient_Calving_PulseGaussianPerturb_', num2str(gauss_mag*10)];
diffu_gauss_title = ['Transient_Calving_DiffuGaussianPerturb_', num2str(gauss_mag*10)];
pulse_gauss_mu_title = ['Transient_Calving_MassUnloading_PulseGaussianPerturb_', num2str(gauss_mag*10)];
diffu_gauss_mu_title = ['Transient_Calving_MassUnloading_DiffuGaussianPerturb_', num2str(gauss_mag*10)];
diffu_gauss_mu_plastic_title = ['Transient_Calving_MassUnloading_DiffuGaussianPerturb_Plastic_', num2str(gauss_mag*10)];

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

% model iteration
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
    foldername = ['wavybed_models_' modelname '/model' identifier];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end

    % step iteration
    for steps = 1:5

        % Cluster parameters
        cluster = generic('name', oshostname(), 'np', 5);
        cluster.interactive = 1;
        waitonlock = 10;

        % Run steps
        org=organizer('repository',foldername,'prefix',['MISMIP_' modelname],'steps',steps);

        if perform(org, 'Mesh_and_Parameterization') % {{{1  STEP 1

            md = model();
            md=parameterize(md, ['wavybed_par_files/' par_name '.par']);
            
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

        if perform(org, 'Transient_Calving_MassUnloading')% {{{1 STEP 5 
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
            k_budd = fric_coef_conversion(law_from, law_to, md, C0, H0, Zb,1);

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
                    Hi = H0 + deltaH;
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0); % still need deltaH to mask out the first several non-perturb years
                    Hi = H0 + deltaH;
                end
                ocean_mask = results(end).MaskOceanLevelset;
                C = mass_unloading(md, Hi, H0, k_budd, C0, C, ocean_mask, 1);
                % append time and assign
                current_time = md.timestepping.start_time;
                C_add_time = [C; current_time + dt_mu];
                md.friction.C = [md.friction.C, C_add_time];

                % solve
                md = solve(md,'tr');

                % save the new result to a separate var
                new_results = [new_results,md.results.TransientSolution(1)];

                % crop out any step outside step 10 to limit some weird
                % over-stepping from ISSM when re-starting the run
                result_tbl = struct2table(md.results.TransientSolution);
                result_tbl = result_tbl(result_tbl.step <= 10, :);
                md.results.TransientSolution = table2struct(result_tbl);
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

            % Calving
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
            
            % Mass unloading
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
            k_budd = fric_coef_conversion(law_from, law_to, md, C0, H0, Zb,1);

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
                    Hi = H0 + deltaH;
                    C = C0;
                else
                    deltaH = mu_time_mask_interp(it)*(results(end).Thickness - H0); % still need deltaH to mask out the first several non-perturb years
                    Hi = H0 + deltaH;
                end
                ocean_mask = results(end).MaskOceanLevelset;
                C = mass_unloading(md, Hi, H0, k_budd, C0, C, ocean_mask, 1);
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
    end %} testbed iteration
end