%% Histogram of the timeseries correlation distribution
% This study generates a suite of density distribution of the correlation
% between h(t) and terminus/grounding line. This is a continuation of the
% map view study (map_ht_gl_corr.m)

%% Main script
% parameters
md_type = 'expt'; % options: "expt", "ctrl"
geom_type = 'deep'; % options: "deep", "shallow"
ds = 500; % meshgrid spacing

% input
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
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

%% Process the data
% iterate over Deep GL models
switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

% pre-allocate
n_simu = size(folder_dir_groups{geom_i}, 1);
gl_hist = cell(1,n_simu);
front_hist = cell(1,n_simu);
dist_to_front = cell(1,n_simu);

% start extracting data
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    switch md_type
        case 'ctrl'
            load([group.folder{j},'/', group.name{j}, '/', ctrl_name])
            modelname = md.miscellaneous.name;
            load(['analyzed_data/calve_only/ht_calve_model_',modelname(9:end),'.mat'])
        case 'expt' % here we use effective pressure experiment by default
            load([group.folder{j},'/', group.name{j}, '/', expt_name])
            modelname = md.miscellaneous.name;
            load(['analyzed_data/mu_calve/ht_mu_calve_model_',modelname(9:end),'.mat'])
        otherwise 
            warning('unsupported input')
    end
    [W, GL, FC] = parse_modelname(modelname);
    
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
    for k = 1:md.mesh.numberofvertices
        coefs = corrcoef(ht_data.gl', h_mat(:,k));
        gl_corr(k) = coefs(1,2);
        coefs = corrcoef(ht_data.front', h_mat(:,k));
        front_corr(k) = coefs(1,2);
    end
    % crop the data: keep the corr.coef. at grounded ice; sample the
    % middle fast-flowing section
    [gl_corr_grid, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, gl_corr, ds);
    [front_corr_grid, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, front_corr, ds);
    % get the grounding line position at the last timestep
    ocean_mask = md.results.TransientSolution(end).MaskOceanLevelset;
    gl_dist = locate_groundingline(md, ocean_mask); 
    % crop
    wid_factor = 0.3;
    x_crop = x < gl_dist;
    y_crop = y < max(y)/2+W*wid_factor & y > max(y)/2-W*wid_factor;
    gl_corr_grid = gl_corr_grid(y_crop, x_crop);
    front_corr_grid = front_corr_grid(y_crop, x_crop);

    % get distance to ice front for each grid point
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
    dist_to_front{j} = abs(front_x_interp - X_crop);
    
    
    % save to cells
    gl_hist{j} = gl_corr_grid;
    front_hist{j} = front_corr_grid;


    disp([modelname,' is processed.'])

end

%% Making the plot
figure('Position',[100,100,800,800])
tiledlayout(3,3,'TileSpacing','none')
for j = 1:n_simu
    nexttile
    scatter(dist_to_front{j}, front_hist{j},...
        'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on;
    scatter(dist_to_front{j}, gl_hist{j},...
        'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
%     if j == 1 % add legend
%         legend
%     end
    ylim([0.5,1])
    if ismember(j, [1,4,7])
        set(gca,'ytick',[0.5,0.7,0.9]);
    else
        set(gca,'ytick',[]);
    end
    if ismember(j, [7,8,9])
        set(gca,'xtick',[1e4,3e4,5e4]);
    else
        set(gca,'xtick',[]);
    end
end
ax = nexttile(4); ax.YLabel.String = 'Correlation';
ax = nexttile(8); ax.XLabel.String = 'Distance to calving front (m)';

exportgraphics(gcf,['plots/correlation_',md_type,'_',geom_type,'.pdf'],"ContentType","vector")
