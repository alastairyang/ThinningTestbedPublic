%% Check extended, but it is one long simulation
% Parameters
pulse_type = 'Diffu'; % options: "Diffu","Pulse"
geom_type = 'deep'; % options: "deep", "shallow"
L = 60e3; % domain length in meter
ds = 200;
dt = 0.05;
idle_t = 5;
last_t = 21; % final year of the perturbation (including 5 idle years)
nt_relax = 100; % twenty time slices in the extended relaxation run

%% Process data
% simulation groups
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

% shallow fully grounded glacier OR deep GL with floating termini
switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

% select the geometry group of our interest
group = folder_dir_groups{geom_i};

% model names
expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8_Extended.mat'];
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
ctrl_extend_name = 'MISMIP_yangTransient_MassUnloading_Extended.mat';

% iterate over the models of our interest
extended_idx = [1,3,7,9];
n_simu = size(folder_dir_groups{geom_i}, 1); % number of simulation in each group
for j = extended_idx

    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_ctrl_extend = load([group.folder{j},'/', group.name{j}, '/', ctrl_extend_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    
    % first, just check the continuity of dH and fric. Coef in the experiment
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    results_tbl_ctrl_extend = struct2table(md_ctrl_extend.results.TransientSolution);
    expt_time = results_tbl_expt.time - results_tbl_expt.time(1);
    ctrl_time = [results_tbl_ctrl.time; results_tbl_ctrl_extend.time] - results_tbl_ctrl.time(1);
    [expt_grid, x, y] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, results_tbl_expt.Thickness, ds);
    [ctrl_grid, ~, ~] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, results_tbl_ctrl.Thickness, ds);
    [ctrl_grid_extend, ~, ~] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, results_tbl_ctrl_extend.Thickness, ds);
    ctrl_grid = cat(1, ctrl_grid, ctrl_grid_extend);
    % get center flow line and normalize
    expt_cfl = transpose(squeeze(expt_grid(:,floor(length(y)/2),:)));
    expt_cfl = expt_cfl - expt_cfl(:,1);
    ctrl_cfl = transpose(squeeze(ctrl_grid(:,floor(length(y)/2),:)));
    ctrl_cfl = ctrl_cfl - ctrl_cfl(:,1);
    % find the shorter time axis between ctrl and expt
    [tmin, tidx] = min([expt_time(end) ctrl_time(end)]);
    if tidx == 1; disp('    Experimental run is shorter'); else; disp('    Control run is shorter');end
    % interpolate onto an evenly spaced time axis
    time = 0:dt:tmin;
    expt_cfl = transpose(interp1(expt_time, expt_cfl', time));
    ctrl_cfl = transpose(interp1(ctrl_time, ctrl_cfl', time));
    dH_cfl   = ctrl_cfl - expt_cfl;
    % expt_dH_cfl = (expt_cfl(:,2:end) - expt_cfl(:,1:end-1))/dt; % dh/dt

    % retrieve grounding line and calving front
    disp('  ........Retrieve calving front and grounding line locations........')
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    t_ctrl_extend = results_tbl_ctrl_extend.time;
    gls_expt = zeros(size(t_expt));
    gls_ctrl = zeros(size(t_ctrl));
    gls_ctrl_extend = zeros(size(t_ctrl_extend));
    front = zeros(size(t_ctrl));
    front_extend = zeros(size(t_ctrl_extend));
    clear results_tbl_ctrl_extend results_tbl_expt results_tbl_ctrl

    for i = 1:length(t_expt)
        % grounding line
        gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);        
    end
    for i = 1:length(t_ctrl)
        gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
        front(i) = locate_calvingfront(md_ctrl, md_ctrl.results.TransientSolution(i).MaskIceLevelset);
    end
    for i = 1:length(t_ctrl_extend)
        gls_ctrl_extend(i) = locate_groundingline(md_ctrl_extend, md_ctrl_extend.results.TransientSolution(i).MaskOceanLevelset);
        front_extend(i) = locate_calvingfront(md_ctrl_extend, md_ctrl_extend.results.TransientSolution(i).MaskIceLevelset);
    end
    disp('  ........Successful!........')
    % interpolate
    gls_ctrl = [gls_ctrl; gls_ctrl_extend];
    front = [front; front_extend];
    t_ctrl = [t_ctrl; t_ctrl_extend];
    gls_ctrl_interp = interp1(t_ctrl-t_ctrl(1), gls_ctrl, time,'spline');
    gls_expt_interp = interp1(t_expt-t_expt(1), gls_expt, time,'spline');
    gls_diff = gls_expt_interp - gls_ctrl_interp;
    front_interp = interp1(t_ctrl-t_ctrl(1), front, time, 'spline');

    % downsample the relaxation period since we do not want to display the full
    % relaxation run
    t_relax_min = last_t+idle_t+dt;
    t_relax_max = max(time);
    t_relax = linspace(t_relax_min, t_relax_max, nt_relax);
    time_new = [time(time<t_relax_min), t_relax];
    % interpolate Delta H, grounding line, and front location to the new
    % relaxation time axis. "dwsp" stands for downsampling
    gls_ctrl_dwsp = interp1(time, gls_ctrl_interp, time_new,'spline');
    gls_diff_dwsp = interp1(time, gls_diff, time_new, 'spline');
    front_dwsp = interp1(time, front_interp, time_new, 'spline');
    dH_cfl_dwsp = transpose(interp1(time', dH_cfl', time_new','spline'));

    % crop out the first 5 years (idle)
    start_t = 5;
    t_keep = time_new>start_t;
    dH_cfl_dwsp = dH_cfl_dwsp(:,t_keep);
    gls_ctrl_dwsp = gls_ctrl_dwsp(t_keep);
    gls_diff_dwsp = gls_diff_dwsp(t_keep);
    front_dwsp = front_dwsp(t_keep);
    time_new = time_new(t_keep);
    time_new = time_new - time_new(1);

    % crop out the fjord (ice is absent)
    plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
    fake_t = 1:length(time_new); % a linear axis for imagesc time axis
    dH_cfl_dwsp = dH_cfl_dwsp(x<=runme_params.terminus0_x,:);
    
    figure('Position',[100,100,700,600]);
    imagesc(fake_t, plot_x, dH_cfl_dwsp,'AlphaData',~isnan(dH_cfl_dwsp));hold on
    colorbar;clim([-5,5])
    yticks_val = [16.5, 26.5, 36.5, 46.5]; set(gca,'YTick',yticks_val);
    set(gca,'YTickLabel',string(runme_params.terminus0_x/1e3 - yticks_val))
    ylabel('Distance to front (km)','FontName','Aria','FontSize',18)
    % add grounding line and calving front
    % unlike imagesc, 'plot' will respect the x value. To make it
    % compatible with imagesc plot, we need a fake x-axis for 'plot'
    current_xaxis = get(gca,'XTick');
    x_max = max(current_xaxis);
    plot(fake_t, gls_ctrl_dwsp/1000,'-k','LineWidth',1.2); hold on
    plot(fake_t, front_dwsp/1000, '-.k','LineWidth',1.2); hold on

    ax = gca; ax.FontSize = 18;
    n_tick_relax = 4; % total N tick labels in the relaxion period
    xticks_p1 = 0:2:last_t; xticks_p2 = floor(linspace(last_t+1, max(time_new), n_tick_relax));
    xticks = [xticks_p1, xticks_p2];

    fake_t_ticks = floor(interp1(time_new, fake_t, xticks));
    TickLabels = [string(xticks)];
    set(gca,'XTick',fake_t_ticks); set(gca,'XTickLabel',TickLabels)
    xtickangle(45)
    % add the vertical dashlines
    xline(interp1(time_new, fake_t, last_t-idle_t),':k','LineWidth',1.6); hold on;
    xline(interp1(time_new, fake_t, last_t),':k','LineWidth',2.4); hold off
    
    % add pulse timeseries plot
    switch pulse_type
        case "Diffu"
            clim([-10,10]);
            [~, pulse, pulse_t] = make_localized_forcing_timeseries();
        case "Pulse"
            clim([-10,10]);
            [pulse, ~, pulse_t] = make_localized_forcing_timeseries();
        otherwise
            error('Unknown type!')
    end
    % change to colormap with contour
    colormap(diverg_colormap(50));
    
    % plot the grounding line as an inset
    p = get(gca, 'Position');
    inset_y = 0.2;
    pp = axes('Parent', gcf, 'Position', [p(1) p(4)*(1-inset_y)+p(1) p(3) p(4)*inset_y*0.95]);
    yyaxis left
    anomaly(fake_t, gls_diff_dwsp);
    set(gca,'XTick',fake_t_ticks)
    set(gca,'XTickLabel',TickLabels)
    xtickangle(45)
    hold on;
    set(gca, 'YDir','reverse')
    xlim([0,interp1(time_new, fake_t, max(time_new))])
    ylabel('$\Delta$ \textsf{GL(m)}','FontSize',18,'Interpreter','latex')
    xlabel('Time (yr)','FontSize',16)
    % add pulse
    yyaxis right
    dtt = runme_params.gauss_timestep;
    pulse_t = pulse_t(start_t/dtt+1:end) - pulse_t(start_t/dtt+1); 
    pulse = pulse(start_t/dtt+1:end);
    tidx2 = time_new<last_t;
    pulse = interp1(pulse_t, pulse, time_new(tidx2));
    plot(fake_t(tidx2), pulse,'k','LineWidth',2); 
    switch pulse_type; case "Diffu"; ylim([0,0.2]); case "Pulse"; ylim([0,1]);otherwise; error('Unknown type');end
    set(gca,'Ytick',[])
    ax = gca;
    ax.YAxis(1).Color = 'k'; 
    ax.YAxis(2).Color = 'k';
    ax.FontSize = 18;
    % add the end-of-perturbation dashline
    xline(interp1(time_new, fake_t, last_t-idle_t),':k','LineWidth',1.6); hold on;
    xline(interp1(time_new, fake_t, last_t),':k','LineWidth',2.4); hold off
    
%     % interp C
%     C = md_expt.friction.C(1:end-1,:);
%     C_time = transpose(md_expt.friction.C(end,:));
%     C_time = C_time - C_time(1);
%     C = transpose(interp1(C_time, C', time));
%     % transform into grid and 
%     C_grid = mat2cell(C',ones(1,length(time)));
%     [C_grid, x, y] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, C_grid, ds);
%     % just the center flowline
%     mid_y = floor(size(C_grid, 2)/2);
%     C_grid = transpose(squeeze(C_grid(:,mid_y,:)));
%     
%     dC = C_grid(:, 2:end) - C_grid(:, 1:end-1);

end

% figure; 
% imagesc(time(1:end-1), x, dC)
% colorbar; clim([-5e2,5e2])
% colormap(diverg_colormap(50))

%% Extended simulation (expt are two separate runs)
% this script is to be excluded from final submission
% here we are checking if the extended run matches the previous local basal pert.

ds = 200;
pulse_type = 'Pulse';
runme_params = readtable('runme_param.csv');

mddir = 'long_models_yang/model_W5000_GL400_FC120000/';

% the perturbation period
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
% the extended period:
%   Both control (ctrl) and experiment (expt) should only contain effective
%   pressure change induced thinning, but the experiment one experiences
%   localized basal perturbation during the previous perturbation period.
ctrl_extend_name = 'MISMIP_yangTransient_MassUnloading_Extended.mat';
expt_extend_name = ['MISMIP_yangTransient_LocalPerturb_',pulse_type,'_Extended.mat'];

md_ctrl = load([mddir ctrl_name]).md;
md_expt = load([mddir expt_name]).md;

results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
% get model information
modelname = md_ctrl.miscellaneous.name;
[W, GL, FC] = parse_modelname(modelname); % width, GL depths, and sliding law coefficient
% isolate the delta H from localized basal perturbation
expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
deltaH = expt_H_interp - [results_tbl_ctrl.Thickness{:}];
deltaH_cell = num2cell(deltaH,1);
[md_grid, x, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
% mask out non-ice part
[mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
md_grid(mask_grid>0) = nan;
md_grid = permute(md_grid,[2,3,1]);

%retrieve grounding line and calving front
t_expt = results_tbl_expt.time;
t_ctrl = results_tbl_ctrl.time;
gls_expt = zeros(size(t_expt));
gls_ctrl = zeros(size(t_ctrl));
front = zeros(size(t_ctrl));
for i = 1:260
    % grounding line
    gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
    gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
    % front
    front(i) = locate_calvingfront(md_ctrl, md_ctrl.results.TransientSolution(i).MaskIceLevelset);
end
% interpolate
gls_expt_interp = interp1(t_expt, gls_expt, t_ctrl);
gls_diff = gls_expt_interp - gls_ctrl;

% crop the initial 5 years no-perturbation period, and some extra
% padding beyond the calving front
start_t = 5; dt = 0.1; end_t = 26;
md_grid = md_grid(:,:,start_t/dt+1:end);
md_grid = md_grid(:,x<=runme_params.terminus0_x,:);
mid_y = floor(size(md_grid,1)/2);
md_grid_mids = squeeze(md_grid(mid_y,:,:));
% crop grounding line vector as well
gls_expt_c = gls_expt_interp(start_t/dt+1:end);
front_c = front(start_t/dt+1:end);
gls_diff_c = gls_diff(start_t/dt+1:end);

% data for plotting
plot_t = 0:0.1:size(md_grid_mids, 2)/10-0.1; last_t = plot_t(end);
plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
fake_t_end = 26;
nt = floor((fake_t_end-plot_t(end))*size(md_grid_mids,2)/plot_t(end));
% get extended data
md_ctrl_extend = load([mddir ctrl_extend_name]).md;
md_expt_extend = load([mddir expt_extend_name]).md;
[rel_GL_extend, front_extend, dH_extend, t_extend, ~, ~, gls_expt_extend] = checkNewSS(md_ctrl_extend, md_expt_extend, ds, nt);
% concatenate to the main data array
fake_t_extend = linspace(plot_t(end),fake_t_end,length(t_extend));
fake_t_extend = fake_t_extend(2:end);
plot_t = [plot_t, fake_t_extend];
gls_diff_c = [gls_diff_c; rel_GL_extend(2:end)];
gls_expt_c = [gls_expt_c; gls_expt_extend(2:end)];
front_c = [front_c; front_extend(2:end)];
md_grid_mids = [md_grid_mids, dH_extend(x<=runme_params.terminus0_x,:)];

% Create the figure (one for each glacier; save all in a folder)
figure('Position',[100,100,700,600]);
imagesc(plot_t, plot_x, md_grid_mids, 'AlphaData',~isnan(md_grid_mids)); hold on
%hax = gca; hax.YTickLabel = flipud(hax.YTickLabel);
yticks_val = [16.5, 26.5, 36.5, 46.5]; set(gca,'YTick',yticks_val);
set(gca,'YTickLabel',string(runme_params.terminus0_x/1e3 - yticks_val))
ylabel('Distance to front (km)','FontName','Aria','FontSize',18)
% add grounding line and calving front
plot(plot_t, gls_expt_c/1000,'-k','LineWidth',1.2); hold on
plot(plot_t, front_c/1000, '-.k','LineWidth',1.2); hold on
% add the end-of-perturbation dashline
xline(16,':k','LineWidth',1.2); hold on
xline(last_t,':k','LineWidth',2.4); hold off
ax = gca; ax.FontSize = 18;
xticks_p1 = 0:2:21; xticks_p2 = 22:2:fake_t_end;
xticks = [xticks_p1, xticks_p2];
fake_t_extend_labels = floor(interp1(fake_t_extend, t_extend(2:end), xticks_p2));
TickLabels = [string(xticks_p1) string(fake_t_extend_labels)];
set(gca,'XTick',xticks); set(gca,'XTickLabel',TickLabels)
xtickangle(45)

% add pulse timeseries plot
switch pulse_type
    case "Diffu"
        clim([-10,10]);
        [~, pulse, pulse_t] = make_localized_forcing_timeseries();
    case "Pulse"
        clim([-5,5]);
        [pulse, ~, pulse_t] = make_localized_forcing_timeseries();
    otherwise
        error('Unknown type!')
end
colormap(diverg_colormap(50));

% plot the grounding line as an inset
p = get(gca, 'Position');
inset_y = 0.2;
pp = axes('Parent', gcf, 'Position', [p(1) p(4)*(1-inset_y)+p(1) p(3) p(4)*inset_y*0.87]);
%plot_gl_t = t_ctrl - t_ctrl(1);
yyaxis left
% ...as an anomaly plot
%anomaly(plot_gl_t(start_t/dt+1:end)-plot_gl_t(start_t/dt+1), gls_diff_c)
anomaly(plot_t, gls_diff_c);
set(gca,'XTick',xticks)
set(gca,'XTickLabel',TickLabels)
xtickangle(45)
hold on;
set(gca, 'YDir','reverse')
xlim([0,xticks(end)])
ylabel('GL(m)','FontSize',18)
xlabel('Time (yr)','FontSize',16)
% add pulse
yyaxis right
dtt = 0.01;
pulse_t = pulse_t(start_t/dtt+1:end) - pulse_t(start_t/dtt+1);
pulse = pulse(start_t/dtt+1:end);
plot(pulse_t, pulse,'k','LineWidth',2);
switch pulse_type; case "Diffu"; ylim([0,0.2]); case "Pulse"; ylim([0,1]);otherwise; error('Unknown type');end
set(gca,'Ytick',[])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.FontSize = 18;
% add the end-of-perturbation dashline
xline(16,':k','LineWidth',1.2); hold on
xline(last_t,':k','LineWidth',2.4); hold off



%% Check the continuity of friction coefficient
ds = 200;
pulse_type = 'Pulse';
runme_params = readtable('runme_param.csv');

mddir = 'long_models_yang/model_W5000_GL400_FC120000/';
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
ctrl_extend_name = 'MISMIP_yangTransient_MassUnloading_Extended.mat';
expt_extend_name = ['MISMIP_yangTransient_LocalPerturb_',pulse_type,'_Extended.mat'];
md_ctrl = load([mddir ctrl_name]).md;
md_expt = load([mddir expt_name]).md;
md_ctrl_extend = load([mddir ctrl_extend_name]).md;
md_expt_extend = load([mddir expt_extend_name]).md;

md_ctrl_C = md_ctrl.friction.C(:,md_ctrl.friction.C(end,:) < md_ctrl_extend.friction.C(end,1));
md_expt_C = md_expt.friction.C(:,md_expt.friction.C(end,:) < md_expt_extend.friction.C(end,1));
% concat time
ctrl_time = [md_ctrl_C(end,:) md_ctrl_extend.friction.C(end,:)];
expt_time = [md_expt_C(end,:) md_expt_extend.friction.C(end,:)];
% concat data
md_ctrl_C = [md_ctrl_C, md_ctrl_extend.friction.C];
md_expt_C = [md_expt_C, md_expt_extend.friction.C];

if expt_time(end) < ctrl_time(end)
    md_ctrl_C = transpose(interp1(ctrl_time', md_ctrl_C', expt_time));
    time = expt_time;
else
    md_expt_C = transpose(interp1(expt_time', md_expt_C', ctrl_time));
    time = ctrl_time;
end

deltaC = md_expt_C(1:end-1,:) - md_ctrl_C(1:end-1,:);
deltaC = mat2cell(deltaC', ones(length(time),1));
[deltaC_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaC', ds);
% just the center flowline
mid_y = floor(size(deltaC_grid, 2)/2);
deltaC_grid = squeeze(deltaC_grid(:,mid_y,:));

%%%%
ctrlC = mat2cell(transpose(md_ctrl_C(1:end-1,:)), ones(length(time),1));
exptC = mat2cell(transpose(md_expt_C(1:end-1,:)), ones(length(time),1));
[ctrlC_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, ctrlC', ds);
[exptC_grid, x, y] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, exptC', ds);
mid_y = floor(size(ctrlC_grid,2)/2);
ctrlC_grid = squeeze(ctrlC_grid(:,mid_y,:));
exptC_grid = squeeze(exptC_grid(:,mid_y,:));

%%%%

% plot
imagesc(time - time(1),x,deltaC_grid')
clim_val = 1e3;
clim([-clim_val, clim_val]); colorbar
