%% Plot h(t) of center flowline from localized basal perturbation
%% The whole center line tabulated along the time axis
ds = 50; % regular meshgrid spacing
pulse_type = 'Diffu'; % types: "Diffu","Pulse"
geom_type = 'shallow'; % types: "deep", "shallow"
expt_type = 'mu'; % types: "no_mu", "mu" (without mass-unloading; with mass unloading)

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
switch expt_type
    case 'no_mu'
        ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
        expt_name = ['MISMIP_yangTransient_Calving_',pulse_type,'GaussianPerturb_8.mat'];
    case 'mu'
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
    otherwise
        error('Unknown experiment type!')
end

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

    % retrieve grounding line and calving front
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    gls_expt = zeros(size(t_expt));
    gls_ctrl = zeros(size(t_ctrl));
    front = zeros(size(t_ctrl));
    for i = 1:260
        % gl
        gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
        gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
        % front
        front(i) = locate_calvingfront(md_ctrl, md_ctrl.results.TransientSolution(i).MaskIceLevelset);
    end
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
    % crop grounding line vector as well
    gls_expt_c = gls_expt(start_t/dt+1:end);
    front_c = front(start_t/dt+1:end);

    % MAKING THE FIGURE!
    ff = figure('Position',[100,100,700,600]);
    plot_t = 0:0.1:size(md_grid_mids, 2)/10-0.1;
    plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
    imagesc(plot_t, plot_x, md_grid_mids, 'AlphaData',~isnan(md_grid_mids)); hold on
    hax = gca; hax.YTickLabel = flipud(hax.YTickLabel);
    set(gca,'YTick',[20,30,40,50])
    %set(gca,'XTick',[])
    ylabel('Distance to front (km)','FontName','Aria','FontSize',15)
    % add grounding line and calving front
    plot(plot_t, gls_expt_c/1000,'-k','LineWidth',1.2); hold on
    plot(plot_t, front_c/1000, '-.k','LineWidth',1.2); hold on
    % add the end-of-perturbation dashline
    xline(16,':k','LineWidth',1.2); hold off
    
    switch pulse_type
        case "Diffu"
            clim([-10,10]);
            [~, pulse, pulse_t] = make_localized_forcing_timeseries();
        case "Pulse"
            clim([-5,5]);
            [pulse, ~, pulse_t] = make_localized_forcing_timeseries();
        otherwise
            error('Unknown type!')
    end
    colormap(diverg_colormap(50));
    
    % plot the grounding line as an inset
    p = get(gca, 'Position'); %cb = colorbar; set(cb,'Position',[cb.Position(1)+0.10 cb.Position(2) cb.Position(3) cb.Position(4)]);
    inset_y = 0.2;
    pp = axes('Parent', gcf, 'Position', [p(1) p(4)*(1-inset_y)+p(1) p(3) p(4)*inset_y]);
    plot_gl_t = t_ctrl - t_ctrl(1);
    yyaxis left
    % ...as an anomaly plot
    anomaly(plot_gl_t(start_t/dt+1:end)-plot_gl_t(start_t/dt+1), gls_diff(start_t/dt+1:end))
    hold on;
    set(gca, 'YDir','reverse')
    xlim([0, end_t-5]); ylabel('GL(m)','FontName','Aria','FontSize',12)
    xlabel('Time (yr)','FontName','Aria')
    % add pulse
    yyaxis right
    dtt = 0.01;
    pulse_t = pulse_t(start_t/dtt+1:end) - pulse_t(start_t/dtt+1); 
    pulse = pulse(start_t/dtt+1:end);
    plot(pulse_t, pulse,'k','LineWidth',2); 
    switch pulse_type; case "Diffu"; ylim([0,0.2]); case "Pulse"; ylim([0,1]);otherwise; error('Unknown type');end
    set(gca,'Ytick',[])
    ax = gca;
    ax.YAxis(1).Color = 'k'; 
    ax.FontSize = 13;
    % add the end-of-perturbation dashline
    xline(16,':k','LineWidth',1.2); hold off
   
    % export 
    plot_name = [md_ctrl.miscellaneous.name(9:end),'_',pulse_type,'_',geom_type,'_',expt_type];
    exportgraphics(gcf, ['plots/localized_stacked_ht/',plot_name,'.png'],'BackgroundColor','none','Resolution',300)
    

end

%% create color bars for the figures, and save as vector files
figure;
imagesc(rand(10,10)); clim([-10,10]);colormap(diverg_colormap(50)); colorbar
exportgraphics(gcf,'plots/colorbar_0_10.pdf','ContentType','vector')

figure;
imagesc(rand(10,10));clim([-5,5]);colormap(diverg_colormap(50)); colorbar
exportgraphics(gcf,'plots/colorbar_0_5.pdf','ContentType','vector')


%% OUDATED: sample selected control points upstream and downstream of the perturbation
gauss_xloc = 3.2e4; % location of center of gaussian perturbation in meter
gcp_ds = 2000; % sampling spacing for ground control points
ds = 50; % regular meshgrid spacing
pulse_type = 'Pulse'; % types: "Diffu","Pulse"
geom_type = 'deep'; % types: "deep", "shallow"
expt_type = 'mu'; % types: "no_mu", "mu" (without mass-unloading; with mass unloading)

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
% plot parameter 
ylabel_i = [1,4,7];
xlabel_i = [7,8,9];
Ws_symb = [10,20,30];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

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
switch expt_type
    case 'no_mu'
        ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
        expt_name = ['MISMIP_yangTransient_Calving_',pulse_type,'GaussianPerturb_8.mat'];
    case 'mu'
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
    otherwise
        error('Unknown experiment type!')
end

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
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
    md_grid(mask_grid>0) = nan;
    % permute
    md_grid = permute(md_grid,[2,3,1]);

    % center flow line
    start_yr = 5;
    end_yr = 21;
    dt = 0.1;
    mid_i = floor(size(md_grid,1)/2);
    crop_ti = start_yr/dt:(end_yr/dt-1);
    % select 5 points upstream and 5 points downstream (spacing 
    [~, gauss_xloc_i] = min(abs(x - gauss_xloc));
    sel_xi = gauss_xloc_i-5*(gcp_ds/ds):(gcp_ds/ds):gauss_xloc_i+5*(gcp_ds/ds);
    sel_xi(sel_xi == gauss_xloc_i) = [];
    % crop out the extra time; select the equally spaced control points
    md_grid_mid = squeeze(md_grid(mid_i,sel_xi,:));
    md_grid_mid = md_grid_mid(:,crop_ti);
    % distance to the perturbation
    sel_x = x(sel_xi);
    dist = sel_x - gauss_xloc; % positive: downstream
    abs_dist = abs(dist);
    unique_dist = unique(abs_dist);
    t_axis = linspace(start_yr, end_yr, size(md_grid_mid,2));
    t_axis = t_axis - min(t_axis); % displaced; relative to start of perturbation year

    color_length = length(unique_dist);
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = cool(color_length);
    colororder(colors_p)
    % color: from nearest to the furthest
    color_axis = transpose(linspace(min(abs_dist), max(abs_dist),length(dist))); % 10^0 to 10^4

    % split upstream and downstream
    ht_down = md_grid_mid(dist > 1000,:);  absdist_down = abs_dist(dist > 1000);
    ht_up   = md_grid_mid(dist <= 1000,:); absdist_up   = abs_dist(dist <= 1000);
    % color interpolated
    color_down = interp1(color_axis, colors_p, absdist_down);
    color_up   = interp1(color_axis, colors_p, absdist_up);

    figure('Position',[100,100,500,700]);
    tiledlayout(2,1,'TileSpacing','none')
    nexttile
    legends = [];
    for i = 1:size(ht_down,1)
        plot(t_axis, ht_down(i,:), 'Color', color_down(i,:)); hold on;
        legends = [legends, num2str(absdist_down(i)/1000) + " km away"];
    end
    title('Downstream')
    xlabel('Time (year)'); ylabel('Elevation change (meter)')
    xlim([0,max(t_axis)]); ylim([-5,10])
    legend(legends,'Location','northwest')
    legends = [];
    nexttile
    for i = 1:size(ht_up,1)
        plot(t_axis, ht_up(i,:), 'Color', color_up(i,:)); hold on;
        legends = [legends, num2str(absdist_up(i)/1000) + " km away"];
    end
    title('Upstream')
    xlabel('Time (year)'); ylabel('Elevation change (meter)')
    xlim([0,max(t_axis)]);ylim([-12,0])
    legend(legends,'Location','southwest')
    % save graphics
    plot_name = [md_ctrl.miscellaneous.name(9:end),'_',pulse_type,'_',geom_type,'_',expt_type];
    exportgraphics(gcf, ['plots/localized_gcp_ht/',plot_name,'.png'],'Resolution',300)
    
end

%% WELL THIS IS Unfinished...(still trying to measure the wave velocity)
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