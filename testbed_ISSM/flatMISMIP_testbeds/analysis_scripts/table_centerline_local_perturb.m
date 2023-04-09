%% The maximum dynamic elevation change magnitude and change rate magnitude
% from the localized basal perturbation experiments
%% The whole center line tabulated along the time axis
ds = 50; % regular meshgrid spacing
dt = 0.1; 
pulse_type = 'Pulse'; % types: "Diffu","Pulse"
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
dH_maxs = zeros(1,n_simu);
dHdt_maxs = zeros(1,n_simu);
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
    clear md_ctrl md_expt
    md_grid(mask_grid>0) = nan;
    md_grid = permute(md_grid,[2,3,1]);

    % retrieve grounding line and calving front
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    clear results_tbl_expt results_tbl_ctrl

    % crop the initial 5 years no-perturbation period, and some extra
    % padding beyond the calving front
    start_t = 5; dt = 0.1; end_t = 26;
    md_grid = md_grid(:,:,start_t/dt+1:end);
    md_grid = md_grid(:,x<=runme_params.terminus0_x,:);
    mid_y = floor(size(md_grid,1)/2);
    md_grid_mids = squeeze(md_grid(mid_y,:,:));

    % find the maximum values
    dH_maxs(j) = max(abs(md_grid_mids),[],'all');
    dHdt = (md_grid_mids(:,2:end) - md_grid_mids(:,1:end-1))/dt;
    dHdt_maxs(j) = max(abs(dHdt),[],'all');

end

% export table as csv
md_names = string(group.name);
tbl = table(md_names,dH_maxs', dHdt_maxs');
tbl.Properties.VariableNames = ["Names","max_dH","max_dHdt"];
tbl_name = ['table_localizedAmp_',pulse_type,'_',geom_type,'_',expt_type];
writetable(tbl, ['result_tables/',tbl_name,'.csv'])
