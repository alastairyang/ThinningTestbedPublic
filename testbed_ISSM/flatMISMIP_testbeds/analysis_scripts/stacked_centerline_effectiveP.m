%% Plot h(t) of center flowline from localized basal perturbation
% The whole center line (no sampling) tabulated along the time axis

%% Plot only the difference between the control and experiment 
% (thickness change attributed to the effective pressure dependence)
gauss_xloc = 3.2e4; % location of center of gaussian perturbation in meter
gcp_ds = 2000; % sampling spacing for ground control points
ds = 50; % regular meshgrid spacing
geom_type = 'deep'; % types: "deep", "shallow"

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);
% read the runme parameters
runme_params = readtable('runme_param.csv');
% colormaps
load('plots/colormap/davos.mat'); load('plots/colormap/lajolla.mat')
lajolla = lajolla(80:end,:);

% split the folder_dir into two groups, separated by grounding line depth
folder_dir_groups = cell(1,2);
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(foldernames_tbl,1),1);
    for j = 1:size(foldernames_tbl.name)
        GL_bool(j) = compare_GLvalue(foldernames_tbl.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = foldernames_tbl(find(GL_bool),:); %#ok<FNDSB> 
end

% we divide the dicussions by the grounding line depth
[~, shallowGL_i] = min(GLs);
[~, deeperGL_i]  = max(GLs);

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';


n_simu = size(folder_dir_groups{geom_i}, 1);
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % get model information
    modelname = md_ctrl.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname);
    % isolate the delta H from localized basal perturbation
    expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Surface{:}]', results_tbl_ctrl.time,'linear','extrap'));
    deltaH = expt_H_interp - [results_tbl_ctrl.Surface{:}];
    deltaH_cell = num2cell(deltaH,1);
    [md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
    % mask out non-ice part
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
    md_grid(mask_grid>0) = nan;
    md_grid = permute(md_grid,[2,3,1]);

    % add the grounding line and signal
    % retrieve grounding line position over time
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    gls_expt = zeros(size(t_expt));
    gls_ctrl = zeros(size(t_ctrl));
    for i = 1:260
        gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
        gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
    end
    % if there is zero (usually the last point), we use the previous GL
    % value
    zero_idx = find(gls_expt == 0); gls_expt(zero_idx) = gls_expt(zero_idx-1);
    zero_idx = find(gls_ctrl == 0); gls_ctrl(zero_idx) = gls_ctrl(zero_idx-1);
    % interpolate
    gls_expt_interp = interp1(t_expt, gls_expt, t_ctrl);
    gls_diff = gls_expt_interp - gls_ctrl;

    % crop the initial 5 years no-perturbation period, and some extra
    % padding beyond the calving front
    start_t = 5; dt = 0.1; end_t = 26;
    md_grid = md_grid(:,:,start_t/dt+1:end);
    md_grid = md_grid(:,x<=runme_params.terminus0_x,:);
    mid_y = floor(size(md_grid,1)/2);
    md_grid_mids = squeeze(md_grid(mid_y,:,:));
    % also crop the grounding line position vector
    gls_expt_c = gls_expt_interp(start_t/dt+1:end);
    gls_ctrl_c = gls_ctrl(start_t/dt+1:end);

    ff = figure('Position',[100,100,700,600]);
    plot_t = 0:0.1:size(md_grid_mids, 2)/10-0.1;
    plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
    imagesc(plot_t, plot_x, md_grid_mids, 'AlphaData',~isnan(md_grid_mids));
    axis off
    colormap(davos);
    clim([-250,0]);
    hold on
    plot(plot_t,gls_ctrl_c/1000,'r-.','LineWidth',1.5); hold on;
    plot(plot_t,gls_expt_c/1000,'k-.','LineWidth',1.5); hold off

    plot_name = [md_ctrl.miscellaneous.name(9:end),'_',pulse_type,'_',geom_type];
    exportgraphics(gcf, ['plots/mu_stacked_ht/',plot_name,'.png'],'BackgroundColor','none','Resolution',300)

end

%% Plot full simulations: H(t) from the control and experiment side-by-side
gauss_xloc = 3.2e4; % location of center of gaussian perturbation in meter
gcp_ds = 2000; % sampling spacing for ground control points
ds = 50; % regular meshgrid spacing
geom_type = 'deep'; % types: "deep", "shallow"

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);
% read the runme parameters
runme_params = readtable('runme_param.csv');
% colormaps
load('plots/colormap/davos.mat'); load('plots/colormap/lajolla.mat')
lajolla = lajolla(80:end,:);

% split the folder_dir into two groups, separated by grounding line depth
folder_dir_groups = cell(1,2);
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(foldernames_tbl,1),1);
    for j = 1:size(foldernames_tbl.name)
        GL_bool(j) = compare_GLvalue(foldernames_tbl.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = foldernames_tbl(find(GL_bool),:); %#ok<FNDSB> 
end

% we divide the dicussions by the grounding line depth
[~, shallowGL_i] = min(GLs);
[~, deeperGL_i]  = max(GLs);

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';


n_simu = size(folder_dir_groups{geom_i}, 1);
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % Get thickness data
    modelname = md_ctrl.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname);
    expt_H = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
    ctrl_H = [results_tbl_ctrl.Surface{:}];
    expt_H_cell = num2cell(expt_H,1); 
    ctrl_H_cell = num2cell(ctrl_H,1);
    [expt_H_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, expt_H_cell, ds);
    [ctrl_H_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, ctrl_H_cell, ds);
    % mask out non-ice part
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
    expt_H_grid(mask_grid>0) = nan;
    ctrl_H_grid(mask_grid>0) = nan;
    expt_H_grid = permute(expt_H_grid,[2,3,1]);
    ctrl_H_grid = permute(ctrl_H_grid,[2,3,1]);
    % Thickness difference wrt the initial thickness
    expt_H_grid = expt_H_grid - expt_H_grid(:,:,1);
    ctrl_H_grid = ctrl_H_grid - ctrl_H_grid(:,:,1);

    % add the grounding line and signal
    % first retrieve grounding line position over time
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    gls_expt = zeros(size(t_expt));
    gls_ctrl = zeros(size(t_ctrl));
    for i = 1:260
        gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
        gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
    end
    % if there is zero (usually the last point), we use the previous GL
    % value
    zero_idx = find(gls_expt == 0); gls_expt(zero_idx) = gls_expt(zero_idx-1);
    zero_idx = find(gls_ctrl == 0); gls_ctrl(zero_idx) = gls_ctrl(zero_idx-1);
    % interpolate
    gls_expt_interp = interp1(t_expt, gls_expt, t_ctrl);
    gls_diff = gls_expt_interp - gls_ctrl;

    % crop the initial 5 years no-perturbation period, and some extra
    % padding beyond the calving front
    start_t = 5; dt = 0.1; end_t = 26;
    expt_H_grid = expt_H_grid(:,:,start_t/dt+1:end);
    ctrl_H_grid = ctrl_H_grid(:,:,start_t/dt+1:end);
    expt_H_grid = expt_H_grid(:,x<=runme_params.terminus0_x,:);
    ctrl_H_grid = ctrl_H_grid(:,x<=runme_params.terminus0_x,:);
    mid_y = floor(size(expt_H_grid,1)/2);
    % center flowline stacked overtime
    expt_H_grid_mids = squeeze(expt_H_grid(mid_y,:,:));
    ctrl_H_grid_mids = squeeze(ctrl_H_grid(mid_y,:,:));
    % also crop the grounding line position vector
    gls_expt_c = gls_expt_interp(start_t/dt+1:end);
    gls_ctrl_c = gls_ctrl(start_t/dt+1:end);

    ff = figure('Position',[100,100,700,600]);
    tiledlayout(1,2,'TileSpacing','none')
    plot_t = 0:0.1:size(expt_H_grid_mids, 2)/10-0.1;
    plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
    [plot_T, plot_X] = meshgrid(plot_t, plot_x);
    nexttile
    contourf(plot_t, fliplr(plot_x), ctrl_H_grid_mids); % fliplr(plot_x) can be interpret as -> distance to ice front
    %imagesc(plot_t, plot_x, ctrl_H_grid_mids, 'AlphaData',~isnan(ctrl_H_grid_mids));
    colormap(davos); clim([-250,0]); hold on;
    plot(plot_t,(runme_params.terminus0_x - gls_ctrl_c)/1000,'r-.','LineWidth',1.5); hold on;
    %contour(plot_T, plot_X, ctrl_H_grid_mids, 5);
    set(gca,'ytick',[])
    nexttile
    contourf(plot_t, fliplr(plot_x), expt_H_grid_mids);
    %imagesc(plot_t, fliplr(plot_x), expt_H_grid_mids, 'AlphaData',~isnan(expt_H_grid_mids))
    colormap(davos); clim([-250,0]); hold on;
    %contour(plot_T, plot_X, expt_H_grid_mids, 5);
    set(gca,'ytick',[])
    colorbar
    hold on
    plot(plot_t,(runme_params.terminus0_x - gls_ctrl_c)/1000,'r-.','LineWidth',1.5); hold on;
    plot(plot_t,(runme_params.terminus0_x - gls_expt_c)/1000,'k-.','LineWidth',1.5); hold off

    plot_name = [md_ctrl.miscellaneous.name(9:end),'_',geom_type];
    exportgraphics(gcf, ['plots/mu_stacked_ht/',plot_name,'.png'],'BackgroundColor','none','Resolution',300)

end

%% OUDATED:WELL THIS IS Unfinished...(still trying to measure the wave velocity)
% next idea is to just select a few points near the perturbation (rather
% than throughout the whole domain right now), if we were to really measure
% the velocity...
locs = cell(size(ht_up,1),1);
for i = 1:size(ht_up,1)
    locs{i} = local_perturb_peaktime(ht_up(i,:),0:0.1:16, "up", 6, 2);
end

%% APPENDIX: Functions
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