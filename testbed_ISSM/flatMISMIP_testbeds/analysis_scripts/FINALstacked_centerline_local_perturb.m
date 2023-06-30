%% Plot thickness change at center flowline from localized basal perturbation
% This creates subplots in figure 3 and 4 of the main texts
% Author: Donglai Yang
% Date: June 28, 2023

%% Parameters
ds = 50; % structured meshgrid spacing
pulse_type = 'Pulse'; % options: "Diffu","Pulse"
geom_type = 'deep'; % options: "deep", "shallow"
expt_type = 'mu'; % options: "no_mu", "mu" (without ice overburden pressure feedback; with ~. In our text, we account for the feedback by default)

%% Read in model parameters
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

% shallow fully grounded glacier OR deep GL with floating termini
switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end
% with ice overburden pressure feedback OR without ~
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

%% Loading the models and processing
n_simu = size(folder_dir_groups{geom_i}, 1); % number of simulation in each group
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % get model information
    modelname = md_ctrl.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname); % width, GL depths, and sliding law coefficient
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
        % grounding line
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

    % Create the figure (one for each glacier; save all in a folder)
    ff = figure('Position',[100,100,700,600]);
    plot_t = 0:0.1:size(md_grid_mids, 2)/10-0.1;
    plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
    imagesc(plot_t, plot_x, md_grid_mids, 'AlphaData',~isnan(md_grid_mids)); hold on
    hax = gca; hax.YTickLabel = flipud(hax.YTickLabel);
    set(gca,'YTick',[20,30,40,50])
    ylabel('Distance to front (km)','FontName','Aria','FontSize',18)
    % add grounding line and calving front
    plot(plot_t, gls_expt_c/1000,'-k','LineWidth',1.2); hold on
    plot(plot_t, front_c/1000, '-.k','LineWidth',1.2); hold on
    % add the end-of-perturbation dashline
    xline(16,':k','LineWidth',1.2); hold off
    ax = gca; ax.FontSize = 18;
    
    % add pulse timeseries plot
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
    p = get(gca, 'Position');
    inset_y = 0.2;
    pp = axes('Parent', gcf, 'Position', [p(1) p(4)*(1-inset_y)+p(1) p(3) p(4)*inset_y]);
    plot_gl_t = t_ctrl - t_ctrl(1);
    yyaxis left
    % ...as an anomaly plot
    anomaly(plot_gl_t(start_t/dt+1:end)-plot_gl_t(start_t/dt+1), gls_diff(start_t/dt+1:end))
    hold on;
    set(gca, 'YDir','reverse')
    xlim([0, end_t-5]); ylabel('GL(m)','FontName','Aria','FontSize',18)
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
    ax.FontSize = 18;
    % add the end-of-perturbation dashline
    xline(16,':k','LineWidth',1.2); hold off
   
    % export 
    plot_name = [md_ctrl.miscellaneous.name(9:end),'_',pulse_type,'_',geom_type,'_',expt_type];
    exportgraphics(gcf, ['plots/localized_stacked_ht/',plot_name,'.png'],'BackgroundColor','none','Resolution',300)
end

%% create color bars for the figures, and save as vector files
% if needed.
% -10 m to 10 m
figure;
imagesc(0*rand(10,10)); clim([-10,10]);colormap(diverg_colormap(50)); colorbar
exportgraphics(gcf,'plots/colorbar_m10_p10.png','Resolution',600)

% -5 m to 5 m
figure;
imagesc(0*rand(10,10));clim([-5,5]);colormap(diverg_colormap(50)); colorbar
exportgraphics(gcf,'plots/colorbar_m5_p5.png','Resolution',600)
