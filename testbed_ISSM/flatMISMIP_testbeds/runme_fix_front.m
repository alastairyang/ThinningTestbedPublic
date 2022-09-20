% read in the table
mdvar_combs = readtable('md_var_combinations.csv');


% iterate over each combination
for jj = 2:size(mdvar_combs,1)

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
    foldername = ['models_' modelname '/model' identifier];
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end

    % RUN
    for steps = 4
        % Mesh
        meshsize = 200;

        % Terminus setup
        terminus0_x = 26500;
        magnitude = 1000;

        % Cluster parameters
        cluster = generic('name', oshostname(), 'np', 2);
        cluster.interactive = 1;
        waitonlock = 10;

        % Run steps
        org=organizer('repository',foldername,'prefix',['MISMIP_' modelname],'steps',steps);

        if perform(org, 'Mesh_and_Parameterization') % {{{1  STEP 1

            md = model();
            md=parameterize(md, ['par_files/' par_name '.par']);

            savemodel(org,md);
        end% }}}
        if perform(org, 'Transient_Steadystate') % {{{1  STEP 2
            % Run the model
            md = loadmodel(org, 'Mesh_and_Parameterization');

            md = setflowequation(md,'SSA','all');

            md.transient.requested_outputs={'default','IceVolume'};

            %md.timestepping = timesteppingadaptive();
            md.timestepping.time_step = 0.010; % default is 0.01
            md.timestepping.start_time = 0;
            md.timestepping.final_time=5; % default is 50
            md.settings.output_frequency=10;

            md.transient.ismovingfront = 0;
            distance_x = md.mesh.x - terminus0_x;
            if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
            expcontourlevelzero(md,distance_x,0,'./TEMP.exp');
            levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp')); % different from distance_x, levelset is all >0 since it takes abs()
            delete('./TEMP.exp');
            md.mask.ice_levelset(md.mesh.x >terminus0_x) = +abs(levelset(md.mesh.x >terminus0_x));
            md.mask.ice_levelset(md.mesh.x<=terminus0_x) = -abs(levelset(md.mesh.x<=terminus0_x));

            md.stressbalance.abstol=NaN;
            md.stressbalance.restol=1e-4;
            md.stressbalance.maxiter=50;
            md.verbose=verbose('convergence',false,'solution',true);
            md.cluster=cluster;
            md.settings.waitonlock=waitonlock;
            md = solve(md,'tr');
            if md.settings.waitonlock == 0
                fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
                fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
                return
            end

            savemodel(org,md);
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

            md.settings.output_frequency = 250;
            md.transient.requested_outputs={'default','IceVolume'};

            md.verbose.solution=1;
            md.cluster = cluster;
            md = solve(md,'tr');

            savemodel(org, md);
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
            
            md.cluster = cluster; generic('name', oshostname(), 'np', 12);
            md = solve(md,'tr');

            savemodel(org, md);
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