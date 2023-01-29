%% Explore map view correlation
% We compute the timeseries correlation at each point on the
% glacier, with respect to both the grounding line and the calving front
% timeseries
% 
% We start by looking at just one model, the one that has the most thinning
% from effective pressure

%% load a selected model from the experiment (mass unloading) 
md_type = 'expt';
md_dir = 'long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat';
load(md_dir)
% load front data; variable name is ht_data
load('analyzed_data/mu_calve/ht_mu_calve_model_W11000_GL400_FC120000');

[front_corr_expt, gl_corr_expt] = gl_front_corr(md, ht_data, md_type);

%% load a selected model from the control
md_type = 'ctrl';
md_dir = 'long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_CalvingOnly.mat';
load(md_dir)
% load front data; variable name is ht_data
load('analyzed_data/calve_only/ht_calve_model_W11000_GL400_FC120000');

[front_corr_ctrl, gl_corr_ctrl] = gl_front_corr(md, ht_data, md_type);

%% APPENDIX: functions
function [front_corr, gl_corr] = gl_front_corr(md, ht_data, md_type)
    
    % put thickness change over time into one big matrix
    result_tbl = struct2table(md.results.TransientSolution);
    n_time = size(result_tbl,1);
    h_mat = [];
    for i = 1:n_time
        h_mat = [h_mat; transpose(cell2mat(result_tbl.Thickness(i)))];
    end
    % correlation
    gl_corr = zeros(md.mesh.numberofvertices, 1);
    front_corr = zeros(md.mesh.numberofvertices, 1);
    for j = 1:md.mesh.numberofvertices
        coefs = corrcoef(ht_data.gl', h_mat(:,j));
        gl_corr(j) = coefs(1,2);
        coefs = corrcoef(ht_data.front', h_mat(:,j));
        front_corr(j) = coefs(1,2);
    end
    % crop the data: keep the grounded ice
    [gl_corr_grid, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, gl_corr, 50);
    [front_corr_grid, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, front_corr, 50);
    mask = md.results.TransientSolution(end).MaskOceanLevelset;
    gl_dist = locate_groundingline(md, mask);
    % crop
    x_crop = x(x < gl_dist);
    gl_corr_grid = gl_corr_grid(:, x < gl_dist);
    front_corr_grid = front_corr_grid(:, x < gl_dist);
    
    % plot
    figure;
    tiledlayout(2,1,'TileSpacing','none','Padding','none')
    nexttile % calving front correlation
    imagesc(x_crop, y, front_corr_grid); colormap copper; clim([0.7,1])
    nexttile % grounding line correlation
    imagesc(x_crop, y, gl_corr_grid); colormap copper; clim([0.7,1])
    cb = colorbar;
    cb.Layout.Tile = 'east';
    % save
    plot_name = ['temporal_correlation_',md_type,'.png'];
    exportgraphics(gcf, ['plots/', plot_name],'Resolution',300)
   
end

%% APPENDIX: un-used code
% CODE BELOW IS NOT USED (but saved for reference)
% %% Correlation between h(t), calving front, and grounding line
% n_gcp = 20;
% sample_interval = 1000; % meter; distance between two control points
% 
% % model parameters and plot parameters
% % read in the model parameter table
% md_vars = readtable('md_var_combinations.csv');
% Ws = sort(unique(md_vars.('fjord_width')));
% GLs = sort(unique(md_vars.('delta_groundingline_depth')));
% FCs = sort(unique(md_vars.('background_friccoef')));
% 
% % specify the scatter plot data symbols
% % GL: circle vs dot; FC: color, light to dark; W: size of the symbol
% Ws_symb = [40,80,110];
% GLs_symb = ["square","o"];
% FCs_symb = [166,32,232;232,32,199;232,32,72]/255;
% 
% ctrl_foldername = 'analyzed_data/calve_only';
% expt_foldername = 'analyzed_data/mu_calve';
% ctrl_folder_prefix = 'ht_calve_';
% expt_folder_prefix = 'ht_mu_calve_';
% ctrl_folder_dir = natsortfiles(dir([pwd '/' ctrl_foldername]));
% expt_folder_dir = natsortfiles(dir([pwd '/' expt_foldername]));
% ctrl_folder_dir = struct2table(ctrl_folder_dir);
% expt_folder_dir = struct2table(expt_folder_dir);
% % remove  '.' and '..'
% bools = cellfun(@(s) ~strcmp(s(1),'.'), ctrl_folder_dir.name);
% ctrl_folder_dir = ctrl_folder_dir(bools,:);
% bools = cellfun(@(s) ~strcmp(s(1),'.'), expt_folder_dir.name);
% expt_folder_dir = expt_folder_dir(bools,:);
% 
% ctrl_folder_dir_groups = cell(1,2);
% for i = 1:length(GLs)
%     % skip the irrelevant ones
%     GL_bool = zeros(size(ctrl_folder_dir,1),1);
%     for j = 1:size(ctrl_folder_dir.name)
%         GL_bool(j) = compare_GLvalue(ctrl_folder_dir.name(j), GLs(i));
%     end
%     % save the respective folder items to a cell
%     ctrl_folder_dir_groups{i} = ctrl_folder_dir(find(GL_bool),:); %#ok<FNDSB> 
% 
% end
% % experiment run
% expt_folder_dir_groups = cell(1,2);
% for i = 1:length(GLs)
%     % skip the irrelevant ones
%     GL_bool = zeros(size(expt_folder_dir,1),1);
%     for j = 1:size(expt_folder_dir.name)
%         GL_bool(j) = compare_GLvalue(expt_folder_dir.name(j), GLs(i));
%     end
%     % save the respective folder items to a cell
%     expt_folder_dir_groups{i} = expt_folder_dir(find(GL_bool),:); %#ok<FNDSB> 
% 
% end
% 
% % we divide the dicussions by the grounding line depth
% [~, shallowGL_i] = min(GLs);
% [~, deeperGL_i]  = max(GLs);
% 
% %% correlation
% % shallower grounding line
% n_simu = size(expt_folder_dir_groups{shallowGL_i}, 1);
% % initialize
% A0s_shallow = zeros(1, n_simu);
% decay_lengths_shallow = zeros(1,n_simu);
% 
% figure
% for j = 1:n_simu
%     % load in data from paths
%     ctrl_path = ctrl_folder_dir_groups{shallowGL_i}(j,:);
%     expt_path = expt_folder_dir_groups{shallowGL_i}(j,:);
%     ctrl = load([ctrl_path.folder{1},'/', ctrl_path.name{1}]);
%     expt = load([expt_path.folder{1},'/', expt_path.name{1}]);
%     ctrl_datas = [ctrl.ht_data.gl(1:end-1)',...
%                   ctrl.ht_data.front(1:end-1)',...
%                   ctrl.ht_data.h];
%     expt_datas = [expt.ht_data.gl(1:end-1)',...
%                   expt.ht_data.front(1:end-1)',...
%                   expt.ht_data.h];
%     ctrl_r = corrcoef(ctrl_datas);
%     ctrl_r = ctrl_r(1:2,3:end);
%     expt_r = corrcoef(expt_datas);
%     expt_r = expt_r(1:2,3:end);
%     plot(1:20, expt_r,'-.k');hold on
%     plot(1:20, ctrl_r,'-.b');hold on;
% end
% 
% %% deeper grounding line
% n_simu = size(expt_folder_dir_groups{deeperGL_i}, 1);
% % initialize
% A0s_deeper = zeros(1, n_simu);
% decay_lengths_deeper = zeros(1,n_simu);
% 
% figure
% ctrl_rs = [];
% expt_rs = [];
% for j = 1:n_simu
%     % load in data from paths
%     ctrl_path = ctrl_folder_dir_groups{deeperGL_i}(j,:);
%     expt_path = expt_folder_dir_groups{deeperGL_i}(j,:);
%     ctrl = load([ctrl_path.folder{1},'/', ctrl_path.name{1}]);
%     expt = load([expt_path.folder{1},'/', expt_path.name{1}]);
%     % look for the first point behind the grounding line 
%     ctrl_point_i = find(ctrl.ht_data.gl(end) > ctrl.ht_data.x,1,'first');
%     expt_point_i = find(expt.ht_data.gl(end) > expt.ht_data.x,1,'first');
%     ctrl_datas = [ctrl.ht_data.gl(1:end-1)',ctrl.ht_data.front(1:end-1)',ctrl.ht_data.h(:,ctrl_point_i)];
%     expt_datas = [expt.ht_data.gl(1:end-1)',expt.ht_data.front(1:end-1)',expt.ht_data.h(:,expt_point_i)];
% 
%     % compare to both grounding line and front positions
% 
%     ctrl_r = corrcoef(ctrl_datas);
%     ctrl_rs = [ctrl_rs; ctrl_r(3,1:2)]; % append the coefficient
%     expt_r = corrcoef(expt_datas);
%     expt_rs = [expt_rs; expt_r(3,1:2)];
% 
%     % plot the thinning time series
%     subplot(2,1,1)
%     plot(ctrl.ht_data.t, ctrl.ht_data.h(:,ctrl_point_i),'-b'); hold on
%     plot(expt.ht_data.t, expt.ht_data.h(:,expt_point_i),'-.k'); hold on
%     
%     % grounding line position
%     subplot(2,1,2)
%     plot(ctrl.ht_data.t, ctrl.ht_data.gl(1:end-1),'-b'); hold on
%     plot(expt.ht_data.t, expt.ht_data.gl(1:end-1),'-.k'); hold on
%     
% end
% subplot(2,1,2)
% plot(ctrl.ht_data.t, ctrl.ht_data.front(1:end-1),'-r')
% 
% figure;
% subplot(1,2,1) % correlation with grounding line
% histogram(ctrl_rs(:,1));hold on; histogram(expt_rs(:,1));
% subplot(1,2,2) % correlation with front
% histogram(expt_rs(:,2));hold on; histogram(expt_rs(:,2))
