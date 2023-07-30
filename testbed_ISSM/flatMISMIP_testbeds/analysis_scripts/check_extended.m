% this script is to be excluded from final submission
% just checking if the extended run matches the previous local basal pert.

ds = 200;
pulse_type = 'Pulse';
runme_params = readtable('runme_param.csv');

dir = 'long_models_yang/model_W5000_GL400_FC120000/';
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
ctrl_extend_name = 'MISMIP_yangTransient_MassUnloading_Extended.mat';
expt_extend_name = ['MISMIP_yangTransient_LocalPerturb_',pulse_type,'_Extended.mat'];
md_ctrl = load([dir ctrl_name]).md;
md_expt = load([dir expt_name]).md;

results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
% get model information
modelname = md_ctrl.miscellaneous.name;
[W, GL, FC] = parse_modelname(modelname); % width, GL depths, and sliding law coefficient
% isolate the delta H from localized basal perturbation
expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
deltaH = expt_H_interp - [results_tbl_ctrl.Thickness{:}];
deltaH_cell = num2cell(deltaH,1);
[md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
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
md_ctrl_extend = load([dir ctrl_extend_name]).md;
md_expt_extend = load([dir expt_extend_name]).md;
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

dir = 'long_models_yang/model_W5000_GL400_FC120000/';
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
ctrl_extend_name = 'MISMIP_yangTransient_MassUnloading_Extended.mat';
expt_extend_name = ['MISMIP_yangTransient_LocalPerturb_',pulse_type,'_Extended.mat'];
md_ctrl = load([dir ctrl_name]).md;
md_expt = load([dir expt_name]).md;
md_ctrl_extend = load([dir ctrl_extend_name]).md;
md_expt_extend = load([dir expt_extend_name]).md;

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
