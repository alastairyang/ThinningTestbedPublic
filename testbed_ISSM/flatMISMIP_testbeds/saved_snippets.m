%% get distance to ice front for each grid point
ice_mask = md.results.TransientSolution(end).MaskIceLevelset;
front_xy = isoline(md, ice_mask,'value',0);
front_y_crop = front_xy.y < max(front_xy.y)/2+W*wid_factor &...
    front_xy.y > max(front_xy.y)/2-W*wid_factor;
front_y = front_xy.y(front_y_crop);
front_x = front_xy.x(front_y_crop);
front_x_interp = interp1(front_y, front_x, transpose(y(y_crop)));
front_x_interp = fillmissing(front_x_interp,"nearest");
front_xy_interp = [front_x_interp, transpose(y(y_crop))];
% distance
[X_crop, Y_crop] = meshgrid(x(x_crop),y(y_crop));
dist_to_front = abs(front_x_interp - X_crop);

%% SNR
    % we use uniform distribution between low and high threshold
    noise_lowamp  = 0; % lower bar 0 m uncertainty
    noise_highamp = 0.5; % higher bar 0.5 m uncertainty
    n_samples = 500; % 200 random sampling of noise amplitude
    amp_rand = (noise_highamp - noise_lowamp).*rand(n_samples,size(expt_S_grid_v,2)) + noise_lowamp;
    noise = transpose(sqrt(amp_rand)/3.*transpose(randn(size(md_grid_v,2),n_samples)));
    snr_ht = zeros(size(md_grid_v,1),n_samples);
    for ni = 1:n_samples
        for k = 1:size(md_grid_v,1)
            if sum(isnan(md_grid_v(k,:)))>0 % skip nan
                snr_ht(k,ni) = nan;
            else
                snr_ht(k,ni) = snr(md_grid_v(k,:), noise(:,ni));
            end
        end
    end


%% Seasonal calving step in runme
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

%% Shear margin weakening step
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

            % Calving
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
            
            % Shear margin weakening

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