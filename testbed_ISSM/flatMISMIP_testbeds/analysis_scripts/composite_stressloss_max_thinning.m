%% force balance
% We find the temporal change in force balance structure
%% Main script
geom_type = "deep"; % options: "deep" or "shallow"
ds = 400; % grid size for the regular grid
sampled_ti = 1:10:240; % sampled time index

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
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

% we divide the dicussions by the grounding line depth
[~, shallowGL_i] = min(GLs);
[~, deeperGL_i]  = max(GLs);
n_simu = size(folder_dir_groups{shallowGL_i}, 1);

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

% save the force balance field data
driving_S_all = cell(2, 260, n_simu);
longi_grad_all = cell(2, 260, n_simu);
later_grad_all = cell(2, 260, n_simu);
basal_R_all = cell(2, 260, n_simu);
gl_x_all = cell(2, 260, n_simu);
front_x_all = cell(2, 260, n_simu);
Ws_md  = zeros(1,n_simu);
GLs_md = zeros(1,n_simu);
FCs_md = zeros(1,n_simu);
dH_max_expt = zeros(1,n_simu);
dH_max_ctrl = zeros(1,n_simu);
dH_sum_expt = zeros(1,n_simu);
dH_sum_ctrl = zeros(1,n_simu);
gl_expt = zeros(1,n_simu);
gl_ctrl = zeros(1,n_simu);

% get both the force balance components
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    index = md_expt.mesh.elements;
    % get the maximum thinning along the centerline
    % experiment
    final_icemask = md_expt.results.TransientSolution(end).MaskIceLevelset;
    final_H = md_expt.results.TransientSolution(end).Thickness; final_H(final_icemask>0) = 0;
    first_H = md_expt.results.TransientSolution(1).Thickness;   first_H(final_icemask>0) = 0;
    deltaH = mesh_to_grid(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, first_H-final_H, ds);
    yi_mid = floor(size(deltaH,1)/2);
    dH_sum_expt(j) = sum(deltaH(yi_mid,:))*ds;
    dH_max_expt(j) = max(deltaH(yi_mid,:),[],'all');
    % control
    final_icemask = md_ctrl.results.TransientSolution(end).MaskIceLevelset;
    final_H = md_ctrl.results.TransientSolution(end).Thickness; final_H(final_icemask>0) = 0;
    first_H = md_ctrl.results.TransientSolution(1).Thickness;   first_H(final_icemask>0) = 0;
    deltaH = mesh_to_grid(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, first_H-final_H, ds);
    yi_mid = floor(size(deltaH,1)/2);
    dH_sum_ctrl(j) = sum(deltaH(yi_mid,:))*ds;
    dH_max_ctrl(j) = max(deltaH(yi_mid,:),[],'all');
    % record model info
    [Ws_md(j), GLs_md(j), FCs_md(j)] = parse_modelname(md_expt.miscellaneous.name);

    % get the maximum grounding line retreat dist
    % first experiment
    final_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(1).MaskOceanLevelset);
    gl_expt(j) = abs(final_gl - first_gl);
    % then control
    final_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(1).MaskOceanLevelset);
    gl_ctrl(j) = abs(final_gl - first_gl);

    % get force balance components
    for ti = sampled_ti
        smooth_L = 2000; % smoothing lengthscale
        [driving_S_expt, basal_R_expt, longi_grad_expt, later_grad_expt, x, y, gl_x_expt, front_x_expt] = calc_force_balance(md_expt,ti,smooth_L);
        [driving_S_ctrl, basal_R_ctrl, longi_grad_ctrl, later_grad_ctrl, ~, ~, gl_x_ctrl, front_x_ctrl] = calc_force_balance(md_ctrl,ti,smooth_L);
        % save
        % experiment
        driving_S_all{1,ti,j} = driving_S_expt; 
        longi_grad_all{1,ti,j} = longi_grad_expt;
        later_grad_all{1,ti,j} = later_grad_expt;
        basal_R_all{1,ti,j} = basal_R_expt;
        gl_x_all{1,ti,j} = gl_x_expt;
        front_x_all{1,ti,j} = front_x_expt;
        % control
        driving_S_all{2,ti,j} = driving_S_ctrl; 
        longi_grad_all{2,ti,j} = longi_grad_ctrl;
        later_grad_all{2,ti,j} = later_grad_ctrl;
        basal_R_all{2,ti,j} = basal_R_ctrl;
        gl_x_all{2,ti,j} = gl_x_ctrl;
        front_x_all{2,ti,j} = front_x_ctrl;

    end
    disp(['model ',num2str(j), ' is completed!'])
end
clear md_ctrl md_expt

% Get integrated stress balance change
runme_params = readtable('runme_param.csv');
front_x = runme_params.terminus0_x;
total_Rs = zeros(2,n_simu);

for ri = 1:2 % run index: first expt, then control
    for j = 1:n_simu % simulations
        ti = sampled_ti(1);
        % get a referential resistive stress estimate at the first timestep
        yi_mid = floor(size(basal_R_all{ri,ti,j},1)/2);
        sample_wid = floor(Ws_md(j)*0.6/ds/2);
        yi_central = yi_mid-sample_wid:yi_mid+sample_wid;
        % calculate the resistive stress
        later_grad_central = later_grad_all{ri,ti,j}(yi_central,:);
        later_grad_mid = mean(later_grad_central,1);
        basal_R_central = basal_R_all{ri,ti,j}(yi_mid,:);
        basal_R_mid = mean(basal_R_central,1);

        % first: we find the integrated resistive stress
        % from gl(end time) to the calving front when basal
        % shear stress is still in place
        gl_x = gl_x_all{ri,sampled_ti(end),j};
        %gl_x = 0;
        front_x = front_x_all{ri,sampled_ti(end),j};
        sample_xi = find(gl_x < x & x < front_x);
        % integrate with trapezoid method
        later_int = trapz(later_grad_mid(sample_xi))*ds;
        basal_int = trapz(basal_R_mid(sample_xi))*ds;
        total_R_ref = later_int+basal_int;

        % at the last timestep, get an estimate of loss of total resistive
        % stress
        ti = sampled_ti(end);
        gl_x = gl_x_all{ri,ti,j};
        %gl_x = 0;
        front_x = front_x_all{ri,ti,j};
        sample_xi = find(gl_x < x & x < front_x);
        later_grad_sample = later_grad_mid(sample_xi);
        later_int = trapz(later_grad_sample)*ds;
        total_R_last = later_int;

        % find the difference
        total_Rs(ri,j) = total_R_ref - total_R_last;
        clear total_R_ref total_R_last

    end
end
%exportgraphics(gcf,'plots/taub_taud_frac.png','Resolution',300)

%% Load narrow + low fric and high fric glaciers
% ... get map view of dH and the locations of max dH
% This and the following sections (until stress loss) can be run
% independently of the previous sections
lowK_md  = load('long_models_yang/model_W5000_GL400_FC30000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
highK_md = load('long_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
param_tbl = readtable('runme_param.csv');

% get dH
lowK_dH = lowK_md.results.TransientSolution(end).Thickness - lowK_md.results.TransientSolution(1).Thickness;
highK_dH = highK_md.results.TransientSolution(end).Thickness - highK_md.results.TransientSolution(1).Thickness;
lowK_mask = lowK_md.results.TransientSolution(end).MaskIceLevelset;
highK_mask = highK_md.results.TransientSolution(end).MaskIceLevelset;
lowK_dH(lowK_mask>0) = nan;
highK_dH(highK_mask>0) = nan;
% get grounding line at last timestep 
lowK_gl = isoline(lowK_md, lowK_md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
highK_gl = isoline(highK_md, highK_md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
[lowK_dH, x, y] = mesh_to_grid(lowK_md.mesh.elements, lowK_md.mesh.x, lowK_md.mesh.y, lowK_dH, 50);
[highK_dH,~, ~] = mesh_to_grid(highK_md.mesh.elements, highK_md.mesh.x, highK_md.mesh.y, highK_dH, 50);
lowK_dH = lowK_dH'; highK_dH = highK_dH';
% crop out the x beyond terminus0_x
xi_keep = find(x < param_tbl.terminus0_x);
x_c = x(xi_keep);
lowK_dH = lowK_dH(xi_keep,:);
highK_dH = highK_dH(xi_keep,:);

% find max dH point along the thalweg
mid_yi = floor(size(lowK_dH,2)/2);
[~, maxi_low]  = max(abs(lowK_dH(:,mid_yi))); 
[~, maxi_high] = max(abs(highK_dH(:,mid_yi)));

%% plot left column: map view of dH and max dH locations
% load colormap
load('plots/colormap/nuuk_polar.mat')
figure('Position',[100,100,350,700])
tiledlayout(1,2,"TileSpacing","none")
nexttile
imagesc(y/1e3,x_c/1e3,lowK_dH); hold on;
scatter(lowK_gl(1).y/1e3, lowK_gl(1).x/1e3,10,'k','filled'); hold on;
colormap(nuuk); clim([-400,0])
scatter(y(mid_yi)/1e3, x_c(maxi_low)/1e3, 100,'r','filled')
ylabel('Along flow distance (km)','FontSize',13)

nexttile
imagesc(y/1e3,x_c/1e3,highK_dH); hold on;
scatter(highK_gl(1).y/1e3, highK_gl(1).x/1e3,10,'k','filled'); hold on;
scatter(y(mid_yi)/1e3, x_c(maxi_high)/1e3, 100,'b','filled')
colormap(nuuk); clim([-400,0])
colorbar

% export
exportgraphics(gcf, 'plots/composite_stressloss/dH.png','Resolution',600);
%% Plot mid-column: H(t), GL, and front
% get h(t) at the max dH point
lowK_dH_pt = plot_select_dhdt(lowK_md, x_c(maxi_low), y(mid_yi));
highK_dH_pt = plot_select_dhdt(highK_md, x_c(maxi_high), y(mid_yi));
% get gl
[lowK_gl,t]  = plot_gl_timeseries(lowK_md);
highK_gl = plot_gl_timeseries(highK_md); 
% get front
lowK_front  = plot_front_timeseries(lowK_md);
highK_front = plot_front_timeseries(highK_md);
t = t - t(1);

% plot
figure('Position',[100,100,350,700])
tiledlayout(2,1,"TileSpacing","none")
nexttile
plot(t,lowK_dH_pt,'-r','LineWidth',2);hold on;
plot(t,highK_dH_pt,'-b','LineWidth',2); hold off;
ylabel('meter','FontSize',13)
legend(["H(t) at low K","H(t) at high K"])
xlim([0,26])
nexttile
plot(t,lowK_gl/1e3,':r','LineWidth',2);hold on;
plot(t,highK_gl/1e3,':b','LineWidth',2); hold on;
plot(t,lowK_front/1e3,'-k','LineWidth',2);hold off; % only one front is needed for plotting
legend(["GL at low K","GL at high K","Calving front"])
xlabel('Year','FontSize',13)
ylabel('meter','FontSize',13)
xlim([0,26])
% export
exportgraphics(gcf,'plots/composite_stressloss/gl_front_ht.png','Resolution',600);

%% Plot frontal resistive stress loss vs GL and max dH 
Ws_symb = [40,100,260];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;
figure('Position',[100,100,350,700]);
tiledlayout(2,1,'TileSpacing','compact')

nexttile % GL retreat vs max thinning; stress loss vs max thinning
% only for the three glaciers
lowK_i  = find((Ws_md == 5e3) + (FCs_md == 0.3e5) == 2);
midK_i  = find((Ws_md == 5e3) + (FCs_md == 0.6e5) == 2);
highK_i = find((Ws_md == 5e3) + (FCs_md == 1.2e5) == 2);
% plot stress loss vs dH
yyaxis left
scatter(dH_max_expt(lowK_i), total_Rs(1,lowK_i)/1e6, 200, FCs_symb(1,:), 'filled'); hold on
scatter(dH_max_expt(midK_i), total_Rs(1,midK_i)/1e6, 200, FCs_symb(2,:), 'filled'); hold on
scatter(dH_max_expt(highK_i), total_Rs(1,highK_i)/1e6, 200, FCs_symb(3,:), 'filled'); hold on
ylabel('Frontal sress loss (MPa m)','FontSize',14)
set(gca,'ycolor','k') 
% plot GL vs dH
yyaxis right
scatter(dH_max_expt(lowK_i), gl_expt(lowK_i)/1e3, 200, FCs_symb(1,:),'filled','Marker','diamond'); hold on
scatter(dH_max_expt(midK_i), gl_expt(midK_i)/1e3, 200, FCs_symb(2,:),'filled','Marker','diamond'); hold on
scatter(dH_max_expt(highK_i), gl_expt(highK_i)/1e3, 200, FCs_symb(3,:),'filled','Marker','diamond'); hold on
ylabel('Grounding line retreat (km)','FontSize',14)
xlabel('Maximum thinning (m)','FontSize',14)
set(gca,'ycolor','k') 

nexttile % total thinning vs GL retreat
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(dH_sum_expt(j),gl_expt(j)/1e3, W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(dH_sum_ctrl(j), gl_ctrl(j)/1e3,W_symb,FC_symb,GL_symb);
    hold on
end
xlabel('Total thinning (m^2)')
ylabel('Grounding line retreat (km)')
ax = gca;
ax.FontSize = 14;
%ax.FontName = 'Aria';
exportgraphics(gcf,'plots/composite_stressloss/stressloss.png','Resolution',600);
exportgraphics(gcf,'plots/composite_stressloss/stressloss.pdf','ContentType','vector')

%% OUTDATED: Plotting maximum thinning and frontal resistive stress loss 
Ws_symb = [40,100,260];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;
figure('Position',[100,100,1100,600]);
tiledlayout(1,3,'TileSpacing','compact')

nexttile % total thinning vs GL retreat
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(gl_expt(j)/1e3, dH_sum_expt(j),W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(gl_ctrl(j)/1e3, dH_sum_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
ylabel('Total thinning (m^2)')
xlabel('Grounding line retreat (km)')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Aria';

nexttile % Max thinning vs GL retreat
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(gl_expt(j)/1e3, dH_max_expt(j), W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(gl_ctrl(j)/1e3, dH_max_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
ylabel('Maximum thinning (m)')
xlabel('Grounding line retreat (km)')
set(gca,'XTick',[2,6,10,14]);
ax = gca; ax.FontSize = 14; ax.FontName = 'Aria';

nexttile % Total stress loss vs max thinning
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(total_Rs(1,j)/1e9, dH_max_expt(j), W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(total_Rs(2,j)/1e9, dH_max_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
ylabel('Maximum thinning (m)','FontName','Aria','FontSize',14)
xlabel('Total resistive stress loss (GPa m)','FontName','Aria','FontSize',14)
ax = gca; ax.FontSize = 14; ax.FontName = 'Aria';

exportgraphics(gcf,'plots/two_thinning_proxies.png','Resolution',600)
