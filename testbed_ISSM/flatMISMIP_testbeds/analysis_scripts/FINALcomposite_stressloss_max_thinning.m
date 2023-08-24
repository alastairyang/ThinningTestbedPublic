%% thinning, GL, and frontal resistive stress loss
% This creates figure 7 in the main text
% Author: Donglai Yang
% Date: June 27, 2023
clear;clc;

%% Stress balance estimation
% first, we need to calculate the frontal resistive stress
% This script does it for all glaciers
geom_type = "deep"; % options: "deep" or "shallow"
ds = 100; % grid size for the regular grid
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

% pre-allocate array
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
    % load the model
    group = folder_dir_groups{geom_i};
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    index = md_expt.mesh.elements;
    % get the maximum thinning along the centerline
    % for experiment
    final_icemask = md_expt.results.TransientSolution(end).MaskIceLevelset;
    final_H = md_expt.results.TransientSolution(end).Thickness; final_H(final_icemask>0) = 0;
    first_H = md_expt.results.TransientSolution(1).Thickness;   first_H(final_icemask>0) = 0;
    deltaH = mesh_to_grid(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, first_H-final_H, ds);
    yi_mid = floor(size(deltaH,1)/2);
    dH_sum_expt(j) = sum(deltaH(yi_mid,:))*ds;
    dH_max_expt(j) = max(deltaH(yi_mid,:),[],'all');
    % for control
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
    % first for experiment
    final_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(1).MaskOceanLevelset);
    gl_expt(j) = abs(final_gl - first_gl);
    % then for control
    final_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(1).MaskOceanLevelset);
    gl_ctrl(j) = abs(final_gl - first_gl);

    % get force balance components
    for ti = sampled_ti
        smooth_L = 2000; % smoothing lengthscale
        [driving_S_expt, basal_R_expt, longi_grad_expt, later_grad_expt, x, y, gl_x_expt, front_x_expt] = calc_force_balance(md_expt,ti,smooth_L,ds);
        [driving_S_ctrl, basal_R_ctrl, longi_grad_ctrl, later_grad_ctrl, ~, ~, gl_x_ctrl, front_x_ctrl] = calc_force_balance(md_ctrl,ti,smooth_L,ds);
        % save!
        %  for experiment
        driving_S_all{1,ti,j} = driving_S_expt; 
        longi_grad_all{1,ti,j} = longi_grad_expt;
        later_grad_all{1,ti,j} = later_grad_expt;
        basal_R_all{1,ti,j} = basal_R_expt;
        gl_x_all{1,ti,j} = gl_x_expt;
        front_x_all{1,ti,j} = front_x_expt;
        %  for control
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

% laod simulation parameters
runme_params = readtable('runme_param.csv');
front_x = runme_params.terminus0_x;
% pre-allocate frontal resistive stress array
frontal_Rs = zeros(2,n_simu);
later_ref = zeros(2,n_simu);
later_post = zeros(2,n_simu);
basal_ref = zeros(2,n_simu);
basal_post = zeros(2,n_simu);
longi_later_ratio = zeros(2,n_simu);
dL = zeros(2,n_simu);

% Compare longitudinal resist. stress to lateral ~
lvl_ratio = cellfun(@(x, y) x./y, longi_grad_all, later_grad_all, 'UniformOutput', false);

%% ------------------------------------------------------------------------
% Frontal stress loss
% From the calculated stress components, we calculate the integrated
% ...frontal resistive stress loss. The 'front' is defined as the distance
% from the last position of the grounding line to the initial position of
% the calving front.
for ri = 1:2 % run index: first expt, then control
    for j = 1:n_simu % testbeds
        ti = sampled_ti(1);
        % get a referential resistive stress estimate at the first timestep
        yi_mid = floor(size(basal_R_all{ri,ti,j},1)/2);
        sample_wid = floor(Ws_md(j)*0.2/ds/2);
        yi_central = yi_mid-sample_wid:yi_mid+sample_wid;

        % first time step, however we use the final position of grounding
        % line
        gl_x = gl_x_all{ri,sampled_ti(end),j};
        later_grad_central = later_grad_all{ri,ti,j}(yi_central,:);
        later_grad_mid = mean(later_grad_central,1);
        basal_R_central = basal_R_all{ri,ti,j}(yi_mid,:);
        basal_R_mid = mean(basal_R_central,1);
        longi_grad_central = longi_grad_all{ri,ti,j}(yi_central,:);
        longi_grad_mid = mean(longi_grad_central,1);
        % determine the effective front location via the finding locations
        % of all first NaN in the center line data
        search_start_idx = 100; % the first dozens are also nan; we avoid them
        longi_xidx = find(isnan(longi_grad_mid(search_start_idx:end)),1,'first')+search_start_idx;
        later_xidx = find(isnan(later_grad_mid(search_start_idx:end)),1,'first')+search_start_idx;
        front_x = min([longi_xidx,later_xidx])*ds;
        sample_xi = find(gl_x < x & x < front_x);
        % integrate with trapezoid method
        later_grad_mid_s = later_grad_mid(sample_xi);
        later_grad_mid_s(isnan(later_grad_mid_s)) = 0.0;
        later_int = trapz(later_grad_mid_s)*ds;
        basal_int = trapz(basal_R_mid(sample_xi))*ds;
        % save needed force component
        later_ref(ri,j) = later_int;
        basal_ref(ri,j) = basal_int;
        % longitudinal resistive stress change
        longi_grad_mid_s = longi_grad_mid(sample_xi);
        longi_grad_mid_s(isnan(longi_grad_mid_s)) = 0.0;
        longi_int = -1*trapz(longi_grad_mid_s)*ds;
        total_R_ref = later_int + basal_int + longi_int;

        % same, but at the last time step
        ti = sampled_ti(end);
        gl_x = gl_x_all{ri,ti,j};
        % get their center flow line data first to find new eff front location
        later_grad_central = later_grad_all{ri,ti,j}(yi_central,:);
        later_grad_mid = mean(later_grad_central,1);
        longi_grad_central = longi_grad_all{ri,ti,j}(yi_central,:);
        longi_grad_mid = mean(longi_grad_central,1);
        longi_xidx = find(isnan(longi_grad_mid(search_start_idx:end)),1,'first')+search_start_idx;
        later_xidx = find(isnan(later_grad_mid(search_start_idx:end)),1,'first')+search_start_idx;
        front_x = min([longi_xidx,later_xidx])*ds;
        sample_xi = find(gl_x < x & x < front_x);
        % MY LATEST CORRECTION
        later_grad_sample = later_grad_mid(sample_xi);
        later_int = trapz(later_grad_sample)*ds;
        basal_R_central = basal_R_all{ri,ti,j}(yi_mid,:);
        basal_R_mid = mean(basal_R_central,1);
        basal_int = trapz(basal_R_mid(sample_xi))*ds;
        % save needed force component
        later_post(ri,j) = later_int;
        basal_post(ri,j) = basal_int;
        dL(ri,j) = front_x - gl_x;
        % longitudinal
        longi_grad_mid_s = longi_grad_mid(sample_xi);
        longi_grad_mid_s(isnan(longi_grad_mid_s)) = 0.0;
        longi_int = -1*trapz(longi_grad_mid_s)*ds;

        % no need to add basal resistive stress since it is the frontal
        % region (all zero)
        total_R_last = later_int + longi_int + basal_int;

        % find the difference (stress loss!)
        frontal_Rs(ri,j) = -1*(total_R_last - total_R_ref);

        % longitudinal resist. stress to lateral ~ ratio in the front
        ratio = lvl_ratio{ri,ti,j};
        ratio_cfl = ratio(yi_mid,:);
        longi_later_ratio(ri,j) = nanmean(ratio_cfl(sample_xi));

    end
end
% estimate the mean lateral resistive stress (in the unit of kPa)
mean_later_ref = later_ref./dL/1e3;


%% Marker specs for the plot
Ws_symb = [40,100,260];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;
% create figure
figure('Position',[100,100,900,600]);
t = tiledlayout(2,2,'TileSpacing','loose');

% ---------------  total thinning vs GL retreat & stress loss ----------------------

axis_bot = nexttile([2,1]); 
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(axis_bot, gl_expt(j)/1e3, dH_sum_expt(j)/1e6, W_symb,FC_symb,'Marker','^','LineWidth',2)
    hold on
end
axis_bot.XAxisLocation = 'bottom';
axis_bot.Box = 'off';
xlabel(axis_bot,'Grounding line retreat (km)','FontSize',15)
axis_bot.FontSize = 14;
set(axis_bot,'YTick',3:6)
ylabel(axis_bot,'Total thinning (km^2)','FontSize',15)
% linear regression
xaxis = min(gl_expt/1e3):max(gl_expt/1e3)*0.05:max(gl_expt/1e3);
Htot_coefs = polyfit(gl_expt/1e3, dH_sum_expt/1e6,1);
line_GHtot = Htot_coefs(1)*xaxis+ Htot_coefs(2);
plot(axis_bot, xaxis, line_GHtot, ':r','LineWidth',1.8); hold on;

% Resistive stress
axis_top = axes('Position', axis_bot.Position);
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(axis_top, frontal_Rs(1,j)/1e6, dH_sum_expt(j)/1e6, W_symb,FC_symb,'filled')
    hold on
end
%axis_top.Color = 'none';
xlabel('Frontal resistive stress loss (MPa m)')
axis_top.FontSize = 14;
axis_top.Color = 'none';
axis_top.XAxisLocation = 'top';
set(axis_top,'YTick',3000:2000:9000)
xlabel(axis_top,'Frontal resistive stress loss (MPa m)','FontSize',15)
set(axis_top,'YTick',[])
set(axis_top,'ycolor','k')
% linear regression
xaxis = min(frontal_Rs(1,:)/1e6):max(frontal_Rs(1,:)/1e6)*0.05:max(frontal_Rs(1,:)/1e6);
Htot_coefs = polyfit(frontal_Rs(1,:),dH_sum_expt/1e6,1);
line_RHtot = Htot_coefs(1)*xaxis+ Htot_coefs(2);
%plot(axis_top, xaxis, line_RHtot, '-r','LineWidth',1.2);

% their respective r-squared 
rsq1 = corrcoef(gl_expt'/1e3, dH_sum_expt'/1e6);
rsq2 = corrcoef(frontal_Rs(1,:)'/1e6, dH_sum_expt'/1e6);

% ----------------- max thinning vs GL retreat & stress loss ---------------
axis_bot = nexttile(t);
set(gca,'YTick',[]);
set(gca,'Xtick',[]);
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot grounding line vs max thinning
    scatter(axis_bot,gl_expt(j)/1e3, dH_max_expt(j),W_symb,FC_symb,'^','LineWidth',2);  hold on
end
axis_bot.XAxisLocation = 'bottom';
axis_bot.Box = 'off';
axis_bot.FontSize = 14;
ylabel(axis_bot,'Max thinning (m)','FontSize',15)
xlabel(axis_bot,'Grounding line retreat (km)','FontSize',15)
% add linear regression and the fitted line
xaxis = min(gl_expt/1e3):max(gl_expt/1e3)*0.05:max(gl_expt/1e3);
lm_GH = fitlm(gl_expt/1e3, dH_max_expt);
r_squared_GH = lm_GH.Rsquared.Ordinary;
line_GH = lm_GH.Coefficients.Estimate(2)*xaxis + lm_GH.Coefficients.Estimate(1);
plot(axis_bot, xaxis, line_GH, ':r','LineWidth',1.8); hold on;

% frontal stress loss
axis_top = axes('Position', axis_bot.Position);
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot frontal resistive stress loss vs max thinning
    scatter(axis_top,frontal_Rs(1,j)/1e6, dH_max_expt(j), W_symb,FC_symb,'filled',GL_symb); hold on
end
axis_top.XAxisLocation = 'top';
axis_top.Box = 'off';
axis_top.Color = 'none';
axis_top.FontSize = 14;
set(axis_top,'YTick',3000:2000:9000)
xlabel(axis_top,'Frontal resistive stress loss (MPa m)','FontSize',15)
% linear regression
xaxis = min(frontal_Rs(1,:)/1e6):max(frontal_Rs(1,:)/1e6)*0.05:max(frontal_Rs(1,:)/1e6);
lm_RH = fitlm(frontal_Rs(1,:)/1e6, dH_max_expt);
r_squared_RH = lm_RH.Rsquared.Ordinary;
line_RH = lm_RH.Coefficients.Estimate(2)*xaxis + lm_RH.Coefficients.Estimate(1);
plot(axis_top, xaxis, line_RH, '-r','LineWidth',1.2); hold on;

% their respective r-squared
rsq3 = corrcoef(gl_expt'/1e3, dH_max_expt');
rsq4 = corrcoef(frontal_Rs(1,:)'/1e6, dH_max_expt');

% ----------------------- Zoom into narrow fjord testbeds -----------------
axis_bot = nexttile(t); 
%axis_bot.Parent = tt;
% % draw inset box in the bottom right
% p = get(gca, 'Position');
% pp = axes('Parent', gcf, 'Position', [p(1)+0.4 p(2)+0.04 p(3)*0.3 p(4)*0.3]);
% only for the three narrow fjord glaciers
lowK_i  = find((Ws_md == 5e3) + (FCs_md == 0.3e5) == 2);
midK_i  = find((Ws_md == 5e3) + (FCs_md == 0.6e5) == 2);
highK_i = find((Ws_md == 5e3) + (FCs_md == 1.2e5) == 2);
% plot stress loss vs max dH
scatter(axis_bot,gl_expt(lowK_i)/1e3, dH_max_expt(lowK_i), 200, FCs_symb(1,:),'Marker','^','LineWidth',2,'MarkerFaceColor','none'); hold on
scatter(axis_bot,gl_expt(midK_i)/1e3, dH_max_expt(midK_i), 200, FCs_symb(2,:),'Marker','^','LineWidth',2,'MarkerFaceColor','none'); hold on
scatter(axis_bot,gl_expt(highK_i)/1e3,dH_max_expt(highK_i), 200, FCs_symb(3,:),'Marker','^','LineWidth',2,'MarkerFaceColor','none'); hold on
axis_bot.XAxisLocation = 'bottom';
xlabel(axis_bot,'Grounding line retreat (km)','FontSize',15)
ylabel(axis_bot,'Max thinning (m)','FontSize',15)
axis_bot.Box = 'off';
set(axis_bot,'ycolor','k') 
set(axis_bot,'YTick',150:50:300)
axis_bot.FontSize = 14;
% plot GL vs max dH
axis_top = axes('Position', axis_bot.Position);
scatter(axis_top,frontal_Rs(1,lowK_i)/1e6, dH_max_expt(lowK_i), 200, FCs_symb(1,:), 'filled'); hold on
scatter(axis_top,frontal_Rs(1,midK_i)/1e6, dH_max_expt(midK_i), 200, FCs_symb(2,:), 'filled'); hold on
scatter(axis_top,frontal_Rs(1,highK_i)/1e6,dH_max_expt(highK_i),200, FCs_symb(3,:), 'filled'); hold on
axis_top.FontSize = 14;
set(axis_top,'YTick',3000:2000:9000)
axis_top.Color = 'none';
axis_top.XAxisLocation = 'top';
xlabel(axis_top,'Frontal resistive stress loss (MPa m)','FontSize',15)
set(axis_top,'YTick',[])
set(axis_top,'ycolor','k')


exportgraphics(gcf,'plots/composite_stressloss/stressloss.png','Resolution',600);


