%% Histogram of the timeseries correlation distribution
% This study generates a suite of density distribution of the correlation
% between h(t) and terminus/grounding line. This is a continuation of the
% map view study (map_ht_gl_corr.m)

%% Main script
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

%% Deeper grounding line
% you can specify if you want to plot experiment (e.g., mass unloading) or
% the control run
md_type = 'ctrl'; % options: "expt", "ctrl"
geom_type = 'shallow'; % options: "deep", "shallow"

% iterate over Deep GL models
%tiledlayout(3,3,'TileSpacing','tight')
switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end
figure('Position',[100,100,1000,800]);
n_simu = size(folder_dir_groups{geom_i}, 1);
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    switch md_type
        case 'ctrl'
            load([group.folder{j},'/', group.name{j}, '/', ctrl_name])
            modelname = md.miscellaneous.name;
            load(['analyzed_data/calve_only/ht_calve_model_',modelname(9:end),'.mat'])
        case 'expt' % here we use mass unloading (effective pressure feedback by default
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
    % crop the data: keep the grounded ice
    [gl_corr_grid, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, gl_corr, 50);
    [front_corr_grid, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, front_corr, 50);
    mask = md.results.TransientSolution(end).MaskOceanLevelset;
    gl_dist = locate_groundingline(md, mask);
    % crop
    x_crop = x(x < gl_dist);
    gl_corr_grid = gl_corr_grid(:, x < gl_dist);
    front_corr_grid = front_corr_grid(:, x < gl_dist);

    % make plot
    %tt = nexttile;
    subplot(3,3,j)
    q = histogram(front_corr_grid,'normalization','probability'); hold on;
    h = histogram(gl_corr_grid,'normalization','probability'); hold on
    h.BinEdges = 0.5:0.02:1;
    q.BinEdges = 0.5:0.02:1;
    title(modelname(9:end))
    % add a grounding line + terminus plot as a inset box
    %axes('Position',[.2 .7 .2 .2])
    p = get(gca, 'Position');
    pp = axes('Parent', gcf, 'Position', [p(1)+0.01 p(2)+0.12 p(3)-0.15 p(4)-0.15]);
    yyaxis left
    plot(ht_data.t, ht_data.front(1:end-1)); hold on;
    set(gca,'xtick',[],'ytick',[])
    yyaxis right
    plot(ht_data.t, ht_data.gl(1:end-1));
    set(gca,'xtick',[],'ytick',[])
    %legend(["GL","Front"])

    disp([modelname,' is processed.'])

end
exportgraphics(gcf,['plots/correlation_',md_type,'_',geom_type,'.pdf'],"ContentType","vector")
