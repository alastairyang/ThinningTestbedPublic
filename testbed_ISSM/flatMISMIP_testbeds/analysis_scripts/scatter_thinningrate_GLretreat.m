%% Scatter plot: the relationship between grounding line retreat distance and effective pressure-induced thinning rate
%% Main script
% parameters
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

% pre-allocate
n_simu = size(folder_dir_groups{1}, 1);
gls_expt = zeros(length(GLs), n_simu);
gls_ctrl = zeros(length(GLs), n_simu);
expt_max_deltaH = zeros(length(GLs), n_simu);
ctrl_max_deltaH = zeros(length(GLs), n_simu);
Ws_md  = zeros(length(GLs), n_simu);
FCs_md = zeros(length(GLs), n_simu);
GLs_md = zeros(length(GLs), n_simu);

for q = 1:length(GLs) % shallow, deep
    % start extracting data
    for j = 1:n_simu
        % read the model
        group = folder_dir_groups{q};
        md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
        md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
        results_tbl_expt = struct2table(md_expt.results.TransientSolution);
        results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
        modelname = md_ctrl.miscellaneous.name;

        % Get thickness data
        [Ws(q,j), GLs(q,j), FCs(q,j)] = parse_modelname(modelname);
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
        % get the maximum thinning magnitude along the center flowline
        mid_y = floor(size(expt_H_grid,1)/2);
        expt_max_deltaH(q,j) = abs(min(squeeze(expt_H_grid(mid_y,:,end))));
        ctrl_max_deltaH(q,j) = abs(min(squeeze(ctrl_H_grid(mid_y,:,end))));

        % get total grounding line retreat distance
        gls_expt(q,j) = abs(locate_groundingline(md_expt, md_expt.results.TransientSolution(end).MaskOceanLevelset) -...
                            locate_groundingline(md_expt, md_expt.results.TransientSolution(1).MaskOceanLevelset));
        gls_ctrl(q,j) = abs(locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(end).MaskOceanLevelset) -...
                            locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(1).MaskOceanLevelset));
        
    end
end

%% Make a scatter plot
% symbols
Ws_symb = [40,100,260];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

% shallow glacier
figure('Position',[100,100,600,600])
for q = 1
    for j = 1:n_simu
        W_symb = Ws_symb(Ws_md(q,j)==Ws); % marker size
        GL_symb = GLs_symb(GLs_md(q,j)==GLs); % marker type (square is shallow; circle is deep)
        FC_symb = FCs_symb(FCs_md(q,j)==FCs,:); % color
        % first plot control
        scatter(gls_ctrl(q,j)/1000, ctrl_max_deltaH(q,j), W_symb, FC_symb,GL_symb,'LineWidth',1.5); hold on;
        % then plot experiment
        scatter(gls_expt(q,j)/1000, expt_max_deltaH(q,j), W_symb, FC_symb,'filled',GL_symb); hold on;
    end
end
xlabel('GL retreat distance (km)','FontName','Aria','FontSize',15)
ylabel('Max thinning magnitude (m)','FontName','Aria','FontSize',15)
exportgraphics(gcf,'plots/thinningmag_GLdist_shallow.png','Resolution',400)

% deep glacier
figure('Position',[100,100,600,600])
for q = 2 
    for j = 1:n_simu
        W_symb = Ws_symb(Ws_md(q,j)==Ws); % marker size
        GL_symb = GLs_symb(GLs_md(q,j)==GLs); % marker type (square is shallow; circle is deep)
        FC_symb = FCs_symb(FCs_md(q,j)==FCs,:); % color
        % first plot control
        scatter(gls_ctrl(q,j)/1000, ctrl_max_deltaH(q,j), W_symb, FC_symb,GL_symb,'LineWidth',1.5); hold on;
        % then plot experiment
        scatter(gls_expt(q,j)/1000, expt_max_deltaH(q,j), W_symb, FC_symb,'filled',GL_symb); hold on;
    end
end
xlabel('GL retreat distance (km)','FontName','Aria','FontSize',15)
ylabel('Max thinning magnitude (m)','FontName','Aria','FontSize',15)
exportgraphics(gcf,'plots/thinningmag_GLdist_deep.png','Resolution',400)
