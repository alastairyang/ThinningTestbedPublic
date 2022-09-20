% runme_free_front.m initializes the model with von Mises calving
% parameterization. With different geometries / model configurations, the 
% initialized grounding line positions are different.
instance = 1;

% two matlab instances run simultaneously
full_index = 1:27;
run_instance1 = [1,2,3,10,11,12,19,20,21]; % the ones with the deepest grounding line.
run_instance2 = setdiff(full_index, run_instance1);

% parameters
meshsize = 200;

% Terminus setup
terminus0_x = 29500; % initialized terminus closed to the front
max_stress_grounded = 1e6; % 1 mPa, 1e6
max_stress_floating = 1e5; % 100 kPa, 1e5

% read in the table
mdvar_combs = readtable('md_var_combinations.csv');

% initiate a table where we record run time
tableSize = [size(mdvar_combs,1), 4];
varNames = ["geometry", "runtime(min)", "step","finished time"];
varTypes = ["string", "double", "int8","datetime"];
runtimeTbl = table('Size', tableSize, 'VariableTypes', varTypes, 'VariableNames',varNames);

% time
tic
% iterate over each combination
if instance == 1
    runs = run_instance1;
    tbl_filename = 'runtime_table_1.csv';
else
    runs = run_instance2;
    tbl_filename = 'runtime_table_2.csv';
end
for jj = runs

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
    foldername = ['models_' modelname '_2' '/model' identifier];
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
            md=parameterize(md, ['par_files/' par_name '.par']);
            
            % specify the number of processor based of 1:1000 element ratio
            np = round(md.mesh.numberofelements/1000);
            cluster = generic('name', oshostname(), 'np', np);
            cluster.interactive = 1;

            savemodel(org,md);
        end% }}}

        if perform(org, 'Transient_Steadystate') % {{{1  STEP 2
            % Run the model for a few simulation year with fixed time step.
            % After the model relaxes a little bit, we will extend the
            % initialization for a much longer time.
            md = loadmodel(org, 'Mesh_and_Parameterization');

            md = setflowequation(md,'SSA','all');

            md.transient.requested_outputs={'default','IceVolume'};
            
            % start with fixed timestep
            md.timestepping.time_step = 0.010; % default is 0.01
            md.timestepping.start_time = 0;
            md.timestepping.final_time = 1;
            md.settings.output_frequency = 10;

            distance_x = md.mesh.x - terminus0_x;
            if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
            isoline(md,distance_x,'value',0,'output','./TEMP.exp');
            levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp')); % different from distance_x, levelset is all >0 since it takes abs()
            delete('./TEMP.exp');
            md.mask.ice_levelset(md.mesh.x >terminus0_x) = +abs(levelset(md.mesh.x >terminus0_x));
            md.mask.ice_levelset(md.mesh.x<=terminus0_x) = -abs(levelset(md.mesh.x<=terminus0_x));
        
            % enable the levelset method
            md.transient.ismovingfront = 1;
            md.calving = calvingvonmises();
            md.calving.stress_threshold_groundedice = max_stress_grounded;
            md.calving.stress_threshold_floatingice = max_stress_floating;
            md.calving.min_thickness = 1;
            % fix other fields
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
            md.levelset.spclevelset = NaN*ones(md.mesh.numberofvertices,1);

            % stressbalance solver parameters
            md.stressbalance.abstol=NaN;
            md.stressbalance.restol=1e-4;
            md.stressbalance.maxiter=80;
            md.verbose=verbose('convergence',false,'solution',true);
            
            md.cluster = cluster;
            md.settings.waitonlock=waitonlock;
            md.toolkits.DefaultAnalysis = bcgslbjacobioptions();
            md.settings.solver_residue_threshold = 1e-5; % I set it; Denis didn't. 
            md = solve(md,'tr');
            if md.settings.waitonlock == 0
                fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
                fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
                return
            end

            savemodel(org,md);

            % run time in seconds, print in minutes
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

            md.cluster = cluster; generic('name', oshostname(), 'np', 13);
	        md.cluster.interactive = 1;
            
            start_time = md.timestepping.start_time;
            final_time = md.timestepping.final_time;
            % use adaptive time stepping
            md.timestepping = timesteppingadaptive();
            md.timestepping.start_time = start_time;
            md.timestepping.time_step_min = 0.005;
            md.timestepping.final_time = md.timestepping.start_time + 120;

            md.settings.output_frequency = 100;
            md.transient.requested_outputs={'default','IceVolume'};
            md.stressbalance.requested_outputs={'default','StrainRatexx',...
                                                'StrainRateyy','StrainRatexy',...
                                                'StrainRateeffective', ...
                                                'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};


            md.verbose.solution=1;
            np = round(md.mesh.numberofelements/1000);
            cluster = generic('name', oshostname(), 'np', np);
            cluster.interactive = 1;
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
        end% }}}

        if perform(org, 'Transient_FrontalMeltOnly') % {{{1  STEP 4
            md = loadmodel(org, 'Transient_Steadystate_Extended');
            md = transientrestart(md);

            % seems like some vertices of base are higher than bed
            % we fix it and make sure these are just numeric issues
            % we fix it and make sure these are just numeric issues
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);
            md.geometry.surface(pos) = md.geometry.base(pos) + md.geometry.thickness(pos);

            % enable moving boundary
            md.transient.ismovingfront = 1;

            % forcing
            peak_rate = 80;
            zero_meltingrate1 = [zeros(md.mesh.numberofvertices,1);md.timestepping.start_time];
            zero_meltingrate2 = [zeros(md.mesh.numberofvertices,1);md.timestepping.start_time + 3];
            peak_meltingrate1 = [peak_rate*ones(md.mesh.numberofvertices,1); md.timestepping.start_time + 6];
            peak_meltingrate2 = [peak_rate*ones(md.mesh.numberofvertices,1); md.timestepping.start_time + 10];
            zero_meltingrate3 = [zeros(md.mesh.numberofvertices,1);md.timestepping.start_time + 13];
            md.frontalforcings.meltingrate = [zero_meltingrate1,...
                                              zero_meltingrate2,...
                                              peak_meltingrate1,...
                                              peak_meltingrate2,...
                                              zero_meltingrate3];
            
            % spclevelset: nan. We do not fix the levelset.
            md.levelset.spclevelset = NaN*ones(md.mesh.numberofvertices,1);
            % calving
            md.calving.calvingrate = md.initialization.vel;
            pos = find(md.mask.ocean_levelset < 0);
            md.calving.calvingrate(pos) = 0.0;
%             pos = find(md.geometry.bed > 0);
%             md.calving.calvingrate(pos) = 0.0;
            % add time
            md.timestepping.final_time = md.timestepping.start_time + 2;
            md.settings.output_frequency = 50;
            
            % output
            md.stressbalance.requested_outputs={'default','StrainRatexx',...
                                                'StrainRateyy','StrainRatexy',...
                                                'StrainRateeffective', ...
                                                'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};

            md.cluster = cluster;
            md = solve(md,'tr');

            savemodel(org, md);

            % run time in seconds, print in minutes
            runTime = toc;
            disp(['    Elapsed time is ' num2str(runTime/60) ' minutes, or ' num2str(runTime/3600) ' hours'])
        end

        if perform(org, 'Transient_CalvingOnly')
            md = loadmodel(org, 'Transient_Steadystate_Extended');
            md = transientrestart(md);

            md.timestepping.final_time = md.timestepping.start_time + 20;
            md.settings.output_frequency = 50;
            
            % we fix it and make sure these are just numeric issues
            pos = find(md.mask.ocean_levelset > 0);
            md.geometry.base(pos) = md.geometry.bed(pos);
            md.geometry.surface(pos) = md.geometry.base(pos) + md.geometry.thickness(pos);

            % forcings
            md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);
%             md.calving = calvingvonmises();
%             md.calving.stress_threshold_groundedice = max_stress_grounded;
%             md.calving.stress_threshold_floatingice = max_stress_floating;
%             md.calving.min_thickness = 1;
            md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);

            % allow the front to move
            md.transient.ismovingfront = 1;
            
            % create sequences of terminus position via spclevelset
            levelset0 = md.mask.ice_levelset;
            md.levelset.spclevelset = [];
            md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
            culmu_magnitude = 0;
            
            calving_start = 8;
            calving_end   = 12;
            
            % the year before it starts: still at the original position 
            md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time + calving_start - 1];
            for time = calving_start : calving_end % only the few years in the middle
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
                %signeddistance = move_terminus_levelset(md, levelset, magnitude, +1, true);
                %md.levelset.spclevelset(:,end+1) = [levelset0; time + 1.0];
            end

            
            md.cluster = cluster; generic('name', oshostname(), 'np', 5);
            md = solve(md,'tr');

            savemodel(org, md);

        end
    end
end