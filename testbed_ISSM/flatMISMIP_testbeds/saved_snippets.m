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

%% locate peak ampitude function
function locs = local_perturb_peaktime(data, t, rel_loc, t_crop, period)
%LOCAL_PERTURB_PEAKTIME find the arrival time of the kinematic wave
%initiated by the localized basal perturbation. The we find the time by
%looking for the time where the peak in signal is observed.

    if nargin < 4; t_crop = 7; period = 2; end % chop out first 7 years
    if nargin < 5; period = 2; end

    % chop the initial extra time
    dt = mean(t(2:end) - t(1:end-1));
    n_crop = t_crop/dt;
    data = data(n_crop+1:end);
    t = t(n_crop+1:end);
    % split into multiple periods
    % find the number of complete periods
    n_period = period/dt;
    num_period = floor(length(data)/n_period);
    data = data(1:n_period*num_period);
    t = t(1:n_period*num_period);
    if size(data,1) ~= 1 % reshape into a row vector
        data = data';
    end
    data_periods = reshape(data, n_period, num_period);
    t_period = 0:dt:period-dt;
    % make signal peak positive
    switch rel_loc
        case "up"
            data_periods = data_periods*(-1); % invert so that max is a peak, not trough
        case "down"
            return
    end
    locs = zeros(size(data_periods,2),1);
    for i = 1:size(data_periods, 2)
        [pk, loc] = findpeaks(data_periods(:,i),t_period);
        if length(loc)>1 || isempty(loc)
            disp('Found multiple/zero peaks! Returning...')
            return 
        else
            locs(i) = loc;
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

%% The first two columns in the old stress loss plot
% The first column: map view of thickness chnage and location of max thinning point
% load model and model parameter table
lowK_md  = load('long_models_yang/model_W5000_GL400_FC30000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
highK_md = load('long_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
param_tbl = readtable('runme_param.csv');

% get dH
lowK_dH = lowK_md.results.TransientSolution(end).Thickness - lowK_md.results.TransientSolution(1).Thickness;
highK_dH = highK_md.results.TransientSolution(end).Thickness - highK_md.results.TransientSolution(1).Thickness;
lowK_mask = lowK_md.results.TransientSolution(end).MaskIceLevelset;
highK_mask = highK_md.results.TransientSolution(end).MaskIceLevelset;
lowK_dH(lowK_mask>0) = nan;
highK_dH(highK_mask>0) = nan;
% get grounding line at last timestep 
lowK_gl = isoline(lowK_md, lowK_md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
highK_gl = isoline(highK_md, highK_md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
[lowK_dH, x, y] = mesh_to_grid(lowK_md.mesh.elements, lowK_md.mesh.x, lowK_md.mesh.y, lowK_dH, 50);
[highK_dH,~, ~] = mesh_to_grid(highK_md.mesh.elements, highK_md.mesh.x, highK_md.mesh.y, highK_dH, 50);
lowK_dH = lowK_dH'; highK_dH = highK_dH';
% crop out the x beyond terminus0_x
xi_keep = find(x < param_tbl.terminus0_x);
x_c = x(xi_keep);
lowK_dH = lowK_dH(xi_keep,:);
highK_dH = highK_dH(xi_keep,:);

% find max dH point along the thalweg
mid_yi = floor(size(lowK_dH,2)/2);
[~, maxi_low]  = max(abs(lowK_dH(:,mid_yi))); 
[~, maxi_high] = max(abs(highK_dH(:,mid_yi)));

% plot left column: map view of dH and max dH locations
% load colormap
load('plots/colormap/nuuk_polar.mat')
figure('Position',[100,100,350,700])
tiledlayout(1,2,"TileSpacing","none")
nexttile
imagesc(y/1e3,x_c/1e3,lowK_dH); hold on;
scatter(lowK_gl(1).y/1e3, lowK_gl(1).x/1e3,10,'k','filled'); hold on;
colormap(nuuk); clim([-400,0])
scatter(y(mid_yi)/1e3, x_c(maxi_low)/1e3, 100,'r','filled')
ylabel('Along flow distance (km)','FontSize',15)

nexttile
imagesc(y/1e3,x_c/1e3,highK_dH); hold on;
scatter(highK_gl(1).y/1e3, highK_gl(1).x/1e3,10,'k','filled'); hold on;
scatter(y(mid_yi)/1e3, x_c(maxi_high)/1e3, 100,'b','filled')
colormap(nuuk); clim([-400,0])
set(gca,'YTick',[])
colorbar('north')

% export
exportgraphics(gcf, 'plots/composite_stressloss/dH.png','Resolution',600);

%% The mid column: H(t), GL, and front
% get h(t) at the max dH point
lowK_dH_pt = plot_select_dhdt(lowK_md, x_c(maxi_low), y(mid_yi));
highK_dH_pt = plot_select_dhdt(highK_md, x_c(maxi_high), y(mid_yi));
% get gl
[lowK_gl,t]  = plot_gl_timeseries(lowK_md);
highK_gl = plot_gl_timeseries(highK_md); 
% get front
lowK_front  = plot_front_timeseries(lowK_md);
highK_front = plot_front_timeseries(highK_md);
t = t - t(1);

% plot
figure('Position',[100,100,350,700])
tiledlayout(3,1,"TileSpacing","none")

nexttile
plot(t,lowK_front/1e3,'-k','LineWidth',2);hold off; % only one front is needed for plotting
ylabel('Front location (km)','FontSize',15)
xline(21, ':k','LineWidth',1.2)
xlim([0,26])
set(gca,'YTick',48:4:58)
set(gca,'Xtick',[])
set(gca,'FontSize',14)

nexttile
plot(t,lowK_dH_pt,'-r','LineWidth',2);hold on;
plot(t,highK_dH_pt,'-b','LineWidth',2); hold off;
xline(21, ':k','LineWidth',1.2)
xlim([0,26])
set(gca,'YTick',-250:100:50)
set(gca,'ycolor','k')
set(gca,'Xtick',[])
set(gca,'FontSize',14)

nexttile
plot(t,lowK_gl/1e3,':r','LineWidth',2);hold on;
plot(t,highK_gl/1e3,':b','LineWidth',2); hold on;
ylabel('GL location (km)','FontSize',15)
xline(21, ':k','LineWidth',1.2)
xlim([0,26])
xlabel('Year','FontSize',15)
set(gca,'YTick',40:4:54)
set(gca,'FontSize',14)
% export
exportgraphics(gcf,'plots/composite_stressloss/gl_front_ht.png','Resolution',600);