%% Plot the temporal evolution of localized basal perturbation
% ...for both the transient and diffused pulses

%% Main script
% read the parameter table for runme file
params = readtable('runme_param.csv');
start_t = 0;

pulse_gauss = params.gauss_mag*make_pulse_gauss(params.pulse_gauss_tscale, params.gauss_efold,...
                                         params.gauss_perturb_repeat_tscale, ...
                                         params.gauss_timestep, params.pulse_gauss_tshift);
diffu_gauss = params.gauss_mag*make_diffu_gauss(params.pulse_gauss_tscale, params.diffu_gauss_tscale,...
                                         params.gauss_efold, params.gauss_perturb_repeat_tscale,...
                                         params.gauss_timestep);

% repeat to full perturbation sequences (8 cycles)
total_cycle = params.perturb_duration/params.gauss_perturb_repeat_tscale;
pulse_full = repmat(pulse_gauss, 1, total_cycle);
diffu_full = repmat(diffu_gauss, 1, total_cycle);
% add the no perturb period
no_perturb_t = 0:params.gauss_timestep:params.no_retreat_duration-params.gauss_timestep;
pulse_full = [zeros(size(no_perturb_t)), pulse_full, zeros(size(no_perturb_t))];
diffu_full = [zeros(size(no_perturb_t)), diffu_full, zeros(size(no_perturb_t))];
% actual time axis
gauss_t = 0:params.gauss_timestep:(params.perturb_duration+2*params.no_retreat_duration-params.gauss_timestep);
gauss_t = gauss_t + (start_t + params.gauss_timestep);

% plot
figure('Position',[100,100,1000,200]);
plot(gauss_t, pulse_full,'r','LineWidth',2); hold on
plot(gauss_t, diffu_full,'b','LineWidth',2); hold off
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1])
xlim([0,max(gauss_t)]);ylim([0,1]);
xlabel('Time (yr)','FontName','Aria','FontSize',15);
ylabel('$\alpha$','Interpreter','latex','FontSize',15)
legend(["Transient Pulse","Diffused Pulse"],'FontName','Aria','FontSize',13)
exportgraphics(gcf,'plots/pulse_forcing.pdf','ContentType','vector')

% Plot the two signals separately
% transient pulse
figure('Position',[100,100,1000,200]);
plot(gauss_t, pulse_full,'r','LineWidth',2); hold on
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1])
xlim([0,max(gauss_t)]);ylim([0,1]);
xlabel('Time (yr)','FontName','Aria','FontSize',15);
ylabel('$\alpha$','Interpreter','latex','FontSize',15)
legend("Transient Pulse",'FontName','Aria','FontSize',13)
exportgraphics(gcf,'plots/TP_forcing.pdf','ContentType','vector')

% 
figure('Position',[100,100,1000,200]);
plot(gauss_t, diffu_full,'b','LineWidth',2); hold off
set(gca,'ytick',[0,0.2])
xlim([0,max(gauss_t)]);ylim([0,0.1]);
xlabel('Time (yr)','FontName','Aria','FontSize',15);
ylabel('$\alpha$','Interpreter','latex','FontSize',15)
legend("Diffused Pulse",'FontName','Aria','FontSize',13)
exportgraphics(gcf,'plots/DP_forcing.pdf','ContentType','vector')



