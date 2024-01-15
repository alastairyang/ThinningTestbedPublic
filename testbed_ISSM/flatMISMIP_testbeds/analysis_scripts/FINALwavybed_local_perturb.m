%% Analyze the localized basal perturbation at an underdeepening section of a fractal rough bed
% This creates figure 6 in the main text
% Date: June 28, 2023
% Author: Donglai Yang

%% Parameter
pulse_type = "Pulse";
ds = 100;
% load parameter table
runme_params = readtable('runme_param.csv');

%% load model data
% ctrl: control
% expt: localized basal perturbation
md_ctrl = load('wavybed_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
md_expt = load('wavybed_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_PulseGaussianPerturb_8.mat').md;

% convert to table for ease of data processing
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
md_grid(mask_grid>0) = nan;
md_grid = permute(md_grid,[2,3,1]);

% retrieve grounding line and calving front
t_expt = results_tbl_expt.time;
t_ctrl = results_tbl_ctrl.time;
gls_expt = zeros(size(t_expt));
gls_ctrl = zeros(size(t_ctrl));
front = zeros(size(t_ctrl));
for i = 1:size(md_ctrl.results.TransientSolution,2)
    % grounding line
    gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
    gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
    % front
    front(i) = locate_calvingfront(md_ctrl, md_ctrl.results.TransientSolution(i).MaskIceLevelset);
end
% interpolate onto the same time axis, in case there is a one timestep
% offset in certain simulations
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
gls_expt_c = gls_expt(start_t/dt+1:end);
front_c = front(start_t/dt+1:end);

%% Create the figure
ff = figure('Position',[100,100,700,600]);
plot_t = 0:0.1:size(md_grid_mids, 2)/10-0.1;
plot_x = x(x<=runme_params.terminus0_x)/1000; % in km
% Hovmoller diagram
imagesc(plot_t, plot_x, md_grid_mids, 'AlphaData',~isnan(md_grid_mids)); hold on
hax = gca; hax.YTickLabel = flipud(hax.YTickLabel);
set(gca,'YTick',[20,30,40,50])
ylabel('Distance to front (km)','FontSize',18)
% add grounding line and calving front
plot(plot_t, gls_expt_c/1000,'-k','LineWidth',1.2); hold on
plot(plot_t, front_c/1000, '-.k','LineWidth',1.2); hold on
xticks = 0:2:21;
set(gca,'XTick',xticks)
% add the end-of-perturbation dashline
xline(16,':k','LineWidth',1.2); hold off
ax = gca;
ax.FontSize = 18;

% add the pulse timeseries depending on its type
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
p = get(gca, 'Position'); %cb = colorbar; set(cb,'Position',[cb.Position(1)+0.10 cb.Position(2) cb.Position(3) cb.Position(4)]);
inset_y = 0.2;
pp = axes('Parent', gcf, 'Position', [p(1) p(4)*(1-inset_y)+p(1) p(3) p(4)*inset_y]);
plot_gl_t = t_ctrl - t_ctrl(1);
yyaxis left
% ...as an anomaly plot
anomaly(plot_gl_t(start_t/dt+1:end)-plot_gl_t(start_t/dt+1), gls_diff(start_t/dt+1:end))
ylim([-50,5])
hold on;
set(gca, 'YDir','reverse')
xlim([0, end_t-5]); 
ylabel('$\Delta$ \textsf{GL(m)}','FontSize',18,'Interpreter','latex')
xlabel('Time (yr)','FontName','Aria')
set(gca,'XTick',xticks)


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
xline(16,':k','LineWidth',1.2); hold off

% save
exportgraphics(gcf,'plots/wavybed_local_perturb.png','Resolution',600)
