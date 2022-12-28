%% Analyzing timeseries from inter-annual basal perturbation experiment
perturb_index = 4;
n_gcp = 20;
sample_interval = 1000; % control point interval in meter, consistent with plot_sampled_dhdt
gauss_xloc = 3.2e4; % gaussian patch location from ice divide
glacier_length = 56500; 
perturb_model_name = "MISMIP_yangTransient_Calving_GaussianPerturb_" + num2str(perturb_index) + ".mat";
nopertb_model_name = "MISMIP_yangTransient_CalvingOnly.mat";
perturb_md = load(perturb_model_name).md;
nopertb_md = load(nopertb_model_name).md;

%% Decompose the timeseries
% we also output the positions of the control points
% the positions are the distance to the influx boundary (x=0)
[perturb_ht,~] = plot_sampled_dhdt(perturb_md, n_gcp);
[nopertb_ht,pos] = plot_sampled_dhdt(nopertb_md, n_gcp);
delta_ht = perturb_ht - nopertb_ht;
% crop out the period where there is no inter-annual perturbation
yr = 5; dt = 0.1;
t = yr/dt;
delta_ht = delta_ht(t+1:end-t,:);
STs = detrend(delta_ht);
LTs = delta_ht - STs;

% make a plot of the decomposed trend and inter-annual variability
fig = figure;
colororder(cool(size(LTs,2)))
subplot(1,2,1)
plot(1:size(LTs,1), LTs);
subplot(1,2,2)
plot(1:size(STs,1), STs);

%% find peaks and 
peaks = zeros(n_gcp,1);
for i = 1:n_gcp
    pks = findpeaks(STs(:,i));
    troughs = findpeaks(-1*STs(:,i));
    peaks(i) = mean(pks) + mean(troughs);
end

%% find the control point immmediately downstream of the center of
% the gaussian perturbation patch
perturb_xi = find(pos > gauss_xloc,1,'last') - 3;
peaks_downstream = peaks(1:perturb_xi);
% estimate the decay length scale
options = optimset('PlotFcns',@optimplotfval);
decay_v = [0.00001];
dist_to_perturb = pos(1:perturb_xi) - pos(perturb_xi);
A0 = max(peaks_downstream);
decay_parms = fminsearch(@minimize_attenuation, decay_v, options, A0, dist_to_perturb, peaks_downstream);
% plot the attenuation model and the data
plot_attenuation(decay_parms, A0, dist_to_perturb, peaks_downstream)


%% APPENDIX: functions
function err = minimize_attenuation(v, A, x, data)
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);
    err = sqrt(mean((A_decay' - data).^2));
end

function plot_attenuation(v, A, x, data)
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);

    figure;
    scatter(x, data, 4,'red'); hold on;
    plot(x, A_decay, '-b'); hold off
    legend(["data","model"])

end


