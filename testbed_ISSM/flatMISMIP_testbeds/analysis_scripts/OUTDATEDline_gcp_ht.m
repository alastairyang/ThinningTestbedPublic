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
    

%     % plot all the lines, colored by up/downstream and distance to pertub
%     figure;
%     for i = 1:length(dist)
%         line_color = interp1(color_axis, colors_p, abs_dist(i));
%         if dist(i) >= 0 % downstream
%             line_style = '-';
%         else % upstream
%             line_style = '--';
%         end
%         plot(t_axis, md_grid_mid(i,:),"Color",line_color,'LineStyle',line_style); hold on
%     end

%     % Animate the delta H(t)
%     for i = 1:260
%         imagesc(x,y,squeeze(md_grid(:,:,i)));title(num2str(i/10));hold on;clim([-5,5]);colormap(diverg_colormap(50));colorbar;
%         pause(0.005)
%     end
    
end