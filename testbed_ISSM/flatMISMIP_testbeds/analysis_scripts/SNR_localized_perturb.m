%% Signal-to-Noise study: can we see the dynamic elevation changes (induced
% by localized basal perturbation) from ICESat-2?
% In this study we use signal to Noise ratio to map out regions where
% changes can be seen given a SNR threshold
%% Main script
gauss_xloc = 3.2e4; % location of center of gaussian perturbation in meter
L = 6e4; % flow domain length
gcp_ds = 2000; % sampling spacing for ground control points
ds = 500; % regular meshgrid spacing
qS = 1; % uncertainty magnitude in meter
pulse_type = 'Pulse'; % types: "Diffu","Pulse"
geom_type = 'deep'; % types: "deep", "shallow"
expt_type = 'plastic'; % types: "plastic", "linear","weertman"

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
    case 'plastic'
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading_Plastic.mat';
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_Plastic_8.mat'];
    case 'linear'
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
    case 'weertman'
        ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
        expt_name = ['MISMIP_yangTransient_Calving_',pulse_type,'GaussianPerturb_8.mat'];
    otherwise
        error('Unknown experiment type!')
end

n_simu = size(folder_dir_groups{geom_i}, 1);
% pre-allocate to save the center flowline data
frac_fl = zeros(n_simu, L/ds);

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
    [expt_S_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_expt.Surface, ds);
    [ctrl_S_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.Surface, ds);
    % mask out non-ice part (same for both experiment and control)
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
    expt_S_grid(mask_grid>0) = nan;
    ctrl_S_grid(mask_grid>0) = nan;
    % reshape the 3D array into 2D (1 space + 1 time)
    expt_S_grid = permute(expt_S_grid,[2,3,1]);
    expt_S_grid_v = reshape(expt_S_grid, [size(expt_S_grid,1)*size(expt_S_grid,2), size(expt_S_grid,3)]);
    ctrl_S_grid = permute(ctrl_S_grid,[2,3,1]);
    ctrl_S_grid_v = reshape(ctrl_S_grid, [size(ctrl_S_grid,1)*size(ctrl_S_grid,2), size(ctrl_S_grid,3)]);

    % use the assumed uncertainty amplitude to derive upper and lower bound
    upper_S = ctrl_S_grid_v + qS;
    lower_S = ctrl_S_grid_v - qS;

    out = (expt_S_grid_v > upper_S) + (expt_S_grid_v < lower_S);
    count = sum(out==1,2);
        
    % re-shape back
    count_grid = reshape(count, [size(expt_S_grid,1), size(expt_S_grid,2)]);
    frac_grid  = count_grid/size(ctrl_S_grid_v,2);

    % plot
%     figure('Position',[100,100,1200,300])
%     imagesc(x,y,frac_grid); % use fraction of the duration
%     clim([0,1]); colorbar; 
%     exportgraphics(gcf, ['plots/SNR/', plotname,'.png'], 'Resolution',400)

    plotname = [modelname(9:end),'_',pulse_type,'_',geom_type,'_', expt_type];
    disp(['Model ', plotname,' is complete!'])

    % save center flow line
    yi_mid = floor(length(y)/2);
    frac_fl(j,:) = frac_grid(yi_mid,:);
end

%% Plot the center flow line
plot_md_idx = [1,3,7,9];
line_widths = [2,2,5,5];
light_color = [252,186,3]/255;
darkr_color = [252,119,3]/255;
line_colors = [light_color;darkr_color;light_color;darkr_color];
% Gaussian patch width on the plot
small_halfwid = 0.08*5000*sqrt(5000/11000)*3;
large_halfwid = 0.08*11000*sqrt(11000/11000)*3;
small_x = [gauss_xloc-small_halfwid, gauss_xloc+small_halfwid,...
           gauss_xloc+small_halfwid, gauss_xloc-small_halfwid]/1000;
large_x = [gauss_xloc-large_halfwid, gauss_xloc+large_halfwid,...
           gauss_xloc+large_halfwid, gauss_xloc-large_halfwid]/1000;
box_y = [0,0,1,1];

% linearly interpolate center line data to finer spatial res
ds_fine = 20; % meter
frac_fl_finer = transpose(interp1(x,frac_fl',x(1):ds_fine:x(end)-ds_fine));

figure('Position',[100,100,600,400])
for i = 1:length(plot_md_idx)
    id = plot_md_idx(i);
    plot(x/1000,frac_fl(id,:),'Color',line_colors(i,:),'LineWidth',line_widths(i))
    hold on;
end
% draw the gray band as
box_color = [189,189,189]/255;
h1 = patch(small_x, box_y, box_color); h1.FaceAlpha = 0.5; hold on
h2 = patch(large_x, box_y, box_color); h2.FaceAlpha = 0.3;
xlabel('Along-flow distance (km)','FontName','Aria','FontSize',14)
ylabel('Fraction of simulation duration','FontName','Aria','FontSize',14)
set(gca,'Xtick',[10,30,50]);
set(gca,'Ytick',[0,0.5,1]);
ax=gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
lineplot_name = [pulse_type,'_',geom_type,'_', expt_type,'_qs',num2str(qS)];
exportgraphics(gcf,['plots/signal_visibility/' lineplot_name '.png'],'Resolution',600)

