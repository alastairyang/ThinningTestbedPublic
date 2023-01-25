% This script makes a tiled plot of the thinning magnitude (h at the last
% timestep minus initial h). It shows both delta H in the control run and
% the delta (delta H) between the control and effective pressure feedback
% experiment.

% Warning: this script is quite a mess...when making changes, making sure
% to be clear about control (ctrl), experiment (exp), shallow, and deeper
% grounding line models

glacier_length = 56500; % initial glacier length (x = 0 to calving front)
length_to_keep = 20000; % 20 km behind the calving front to keep for the plots

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));

% specify the scatter plot data symbols
% GL: circle vs dot; FC: color, light to dark; W: size of the symbol
Ws_symb = [40,80,110];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);

plot_idx = 0;

ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

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

% shallower grounding line
n_simu = size(folder_dir_groups{shallowGL_i}, 1);
% pre-allocate
deltaH_ctrl = cell(n_simu, 2);
deltaH_expt = cell(n_simu, 2);
gl_cells_ctrl = cell(n_simu,2);
gl_cells_expt = cell(n_simu,2);

%% Control
% Shallow grounding line
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{shallowGL_i};
    ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    % plot the delta H
    deltaH = ctrl.results.TransientSolution(end).Surface - ...
             ctrl.results.TransientSolution(1).Surface;
    % crop the extents
    % keep: from calving front to 20 km behind
    mask = ctrl.results.TransientSolution(end).MaskIceLevelset;
    front_dist = locate_calvingfront(ctrl, mask);
    upstream_dist = front_dist - length_to_keep;
    [grid_deltaH, x, y] = mesh_to_grid(ctrl.mesh.elements, ctrl.mesh.x, ctrl.mesh.y, deltaH, 50);
    x_crop = x(x < front_dist & x > upstream_dist);
    grid_deltaH_crop = grid_deltaH(:, x < front_dist & x > upstream_dist);

    % save some grounding line positions (sample every other year)
    sampled_t_i = 1:20:size(ctrl.results.TransientSolution,2);
    results_tbl = struct2table(ctrl.results.TransientSolution);
    sampled_t = results_tbl.time(sampled_t_i);
    gl_masks  = results_tbl.MaskOceanLevelset(sampled_t_i);
    % get the x,y of the grounding line from each mask
    xs = [];
    ys = [];
    all_gls = cell(1,length(gl_masks));
    for jj = 1:length(gl_masks)
        all_gls{jj} = isoline(ctrl, gl_masks{jj}, 'value', 0);
    end
    gl_cells_ctrl{j,shallowGL_i}.data = all_gls;
    gl_cells_ctrl{j,shallowGL_i}.t = sampled_t(jj);
    % save to the cells
    deltaH_ctrl{j, shallowGL_i}.data = grid_deltaH_crop;
    deltaH_ctrl{j, shallowGL_i}.x = x_crop;
    deltaH_ctrl{j, shallowGL_i}.y = y;
end

% deeper grounding line
n_simu = size(folder_dir_groups{deeperGL_i}, 1);
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{deeperGL_i};
    ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    % plot the delta H
    deltaH = ctrl.results.TransientSolution(end).Surface - ...
             ctrl.results.TransientSolution(1).Surface;
    % crop the extents
    % keep: from calving front to 20 km behind
    mask = ctrl.results.TransientSolution(end).MaskIceLevelset;
    front_dist = locate_calvingfront(ctrl, mask);
    upstream_dist = front_dist - length_to_keep;
    [grid_deltaH, x, y] = mesh_to_grid(ctrl.mesh.elements, ctrl.mesh.x, ctrl.mesh.y, deltaH, 50);
    x_crop = x(x < front_dist & x > upstream_dist);
    grid_deltaH_crop = grid_deltaH(:, x < front_dist & x > upstream_dist);

    % save some grounding line positions (sample every other year)
    sampled_t_i = 1:20:size(ctrl.results.TransientSolution,2);
    results_tbl = struct2table(ctrl.results.TransientSolution);
    sampled_t = results_tbl.time(sampled_t_i);
    gl_masks  = results_tbl.MaskOceanLevelset(sampled_t_i);
    % get the x,y of the grounding line from each mask
    xs = [];
    ys = [];
    all_gls = cell(1,length(gl_masks));
    for jj = 1:length(gl_masks)
        all_gls{jj} = isoline(ctrl, gl_masks{jj}, 'value', 0);
    end
    gl_cells_ctrl{j,deeperGL_i}.data = all_gls;
    gl_cells_ctrl{j,deeperGL_i}.t = sampled_t(jj);

    % save to the cells
    deltaH_ctrl{j, deeperGL_i}.data = grid_deltaH_crop;
    deltaH_ctrl{j, deeperGL_i}.x = x_crop;
    deltaH_ctrl{j, deeperGL_i}.y = y;
end

%% Experiment (effective pressure feedback)
% Shallow grounding line
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{shallowGL_i};
    expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    % plot the delta H
    deltaH = expt.results.TransientSolution(end).Surface - ...
             expt.results.TransientSolution(1).Surface;
    % crop the extents
    % keep: from calving front to 20 km behind
    mask = expt.results.TransientSolution(end).MaskIceLevelset;
    front_dist = locate_calvingfront(expt, mask);
    upstream_dist = front_dist - length_to_keep;
    [grid_deltaH, x, y] = mesh_to_grid(expt.mesh.elements, expt.mesh.x, expt.mesh.y, deltaH, 50);
    x_crop = x(x < front_dist & x > upstream_dist);
    grid_deltaH_crop = grid_deltaH(:, x < front_dist & x > upstream_dist);

    % save some grounding line positions (sample every other year)
    sampled_t_i = 1:20:size(expt.results.TransientSolution,2);
    results_tbl = struct2table(expt.results.TransientSolution);
    sampled_t = results_tbl.time(sampled_t_i);
    gl_masks  = results_tbl.MaskOceanLevelset(sampled_t_i);
    % get the x,y of the grounding line from each mask
    xs = [];
    ys = [];
    all_gls = cell(1,length(gl_masks));
    for jj = 1:length(gl_masks)
        all_gls{jj} = isoline(expt, gl_masks{jj}, 'value', 0);
    end
    gl_cells_expt{j,shallowGL_i}.data = all_gls;
    gl_cells_expt{j,shallowGL_i}.t = sampled_t(jj);

    % save to the cells
    deltaH_expt{j, shallowGL_i}.data = grid_deltaH_crop;
    deltaH_expt{j, shallowGL_i}.x = x_crop;
    deltaH_expt{j, shallowGL_i}.y = y;
end

% deeper grounding line
n_simu = size(folder_dir_groups{deeperGL_i}, 1);
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{deeperGL_i};
    expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    % plot the delta H
    deltaH = expt.results.TransientSolution(end).Surface - ...
             expt.results.TransientSolution(1).Surface;
    % crop the extents
    % keep: from calving front to 20 km behind
    mask = expt.results.TransientSolution(end).MaskIceLevelset;
    front_dist = locate_calvingfront(expt, mask);
    upstream_dist = front_dist - length_to_keep;
    [grid_deltaH, x, y] = mesh_to_grid(expt.mesh.elements, expt.mesh.x, expt.mesh.y, deltaH, 50);
    x_crop = x(x < front_dist & x > upstream_dist);
    grid_deltaH_crop = grid_deltaH(:, x < front_dist & x > upstream_dist);

    % save some grounding line positions (sample every other year)
    sampled_t_i = 1:20:size(expt.results.TransientSolution,2);
    results_tbl = struct2table(expt.results.TransientSolution);
    sampled_t = results_tbl.time(sampled_t_i);
    gl_masks  = results_tbl.MaskOceanLevelset(sampled_t_i);
    % get the x,y of the grounding line from each mask
    xs = [];
    ys = [];
    all_gls = cell(1,length(gl_masks));
    for jj = 1:length(gl_masks)
        all_gls{jj} = isoline(expt, gl_masks{jj}, 'value', 0);
    end
    gl_cells_expt{j,deeperGL_i}.data = all_gls;
    gl_cells_expt{j,deeperGL_i}.t = sampled_t(jj);

    % save to the cells
    deltaH_expt{j, deeperGL_i}.data = grid_deltaH_crop;
    deltaH_expt{j, deeperGL_i}.x = x_crop;
    deltaH_expt{j, deeperGL_i}.y = y;
end

%% Make tiled plots

figure('Position',[200,200,500,500]);
h = tiledlayout(3,3, 'TileSpacing', 'none', 'Padding', 'none');
for j = 1:n_simu
    nexttile
    imagesc(deltaH_expt{j,shallowGL_i}.x, deltaH_expt{j,shallowGL_i}.y,...
            deltaH_expt{j,shallowGL_i}.data - deltaH_ctrl{j,shallowGL_i}.data);
    colormap("pink")
    clim([-300,0])
    hold on
    % select the data for this glacier
    gl_data = gl_cells_expt{j,shallowGL_i};
    for k = 1:length(gl_data.data)
        scatter(gl_data.data{k}(1).x, gl_data.data{k}(1).y, 5,'filled',... ...
                'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
        hold on
    end
    colororder(cool(length(gl_data.data)))
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end
cb = colorbar;
cb.Layout.Tile = 'east';
exportgraphics(gcf,'plots/dH_from_feedback_shallowGL.pdf')

%%
figure('Position',[800,200,500,500]);
h = tiledlayout(3,3, 'TileSpacing', 'none', 'Padding', 'none');
for j = 1:n_simu
    nexttile
    imagesc(deltaH_expt{j,deeperGL_i}.x, deltaH_expt{j,deeperGL_i}.y,...
            deltaH_expt{j,deeperGL_i}.data - deltaH_ctrl{j,deeperGL_i}.data);
    colormap("pink")
    clim([-350,0])
    hold on;
    % select the data for this glacier
    gl_data = gl_cells_expt{j,deeperGL_i};
    for k = 1:length(gl_data.data)
        scatter(gl_data.data{k}(1).x, gl_data.data{k}(1).y, 5,'filled',... ...
                'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
        hold on
    end
    colororder(cool(length(gl_data.data)))
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])

end
cb = colorbar;
cb.Layout.Tile = 'east';
exportgraphics(gcf,'plots/dH_from_feedback_deepGL.pdf')