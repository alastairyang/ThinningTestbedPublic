%% use diffusion ansatz to fit the data 
% get the dh/dt and x data from simulation
n_gcp = 20;
glacier_length = 56500;
retreat_rate_max = 1000;
perturb_duration = 16;
no_retreat_duration = 5;
[dhdt, pos] = plot_sampled_dhdt(md, n_gcp);
pos = glacier_length - pos;

% retreat sequenece in the simulation
retreat_advance = linspace(100,retreat_rate_max, perturb_duration/2);
retreat_slow = flip(retreat_advance);
retreat_no = zeros(1,no_retreat_duration);
retreat_sequence = [retreat_no, retreat_advance, retreat_slow, retreat_no];
% interp to 0.1 year
old_t = linspace(0,26,length(retreat_sequence));
new_t = linspace(0, 26, size(dhdt,1));
retreat_seq_interp = interp1(old_t, retreat_sequence, new_t);
% culmulative distance
retreat_dist = cumtrapz(new_t, retreat_seq_interp);
%% fit model to data
parms_all = zeros(n_gcp,4);
%parms_all = zeros(n_gcp,3);
for gcp_i = 1:n_gcp
    x = pos(gcp_i);%*ones(size(retreat_dist)) - retreat_dist;
    data = dhdt(:,gcp_i);
    beta_guess = max(data) - min(data);
    disp_guess = min(data);
    v = [1000, 2000000, disp_guess, beta_guess]; % initial guess
    options = optimset('PlotFcns',@optimplotfval);
    parms = fminsearch(@minimize_ansatz, v, options, x, data);
    
    % plot model fit and data
    %plot_ansatz(parms, x, data);hold on;
    parms_all(gcp_i,:) = parms;
end

% Peclet number
dist_terminus = glacier_length - pos;
Pe = parms(:,1).*dist_terminus'./parms_all(:,2);

%% fit a decay model
% we use the multiplier "beta" from the previous model
decay_v = [0.00001];
A0 = parms_all(1,4);
decay_parms = fminsearch(@minimize_attenuation, decay_v, options, A0, pos-pos(1), parms_all(:,4));

% plot
plot_attenuation(decay_parms, A0, pos-pos(1), parms_all(:,4))
%% Ansatz

function err = minimize_ansatz(v, x, data)
    u = v(1); % diffusive wave velocity
    D = v(2); % diffusivity
    b = v(3); % displacement
    beta = v(4); % multiplier

    t = linspace(0,26,length(data)); % 26 years
    x_hat = (x-u*t)./(sqrt(4.*D.*t));
    H = b + beta*erf(x_hat);
    err = sqrt(mean((data - H').^2));
end

function err = minimize_attenuation(v, A, x, data)
    %A = v(1); % initial amplitude
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);
    err = sqrt(mean((A_decay' - data).^2));
end

function plot_ansatz(v, x, data)
    u = v(1); % diffusive wave velocity
    D = v(2); % diffusivity
    b = v(3); % displacement
    beta = v(4); % multiplier

    t = linspace(0,26,length(data)); % 26 years
    x_hat = (x-u*t)./(sqrt(4.*D.*t));
    H = b + beta*erf(x_hat);
    figure;
    plot(t, data,'-.b');hold on
    plot(t, H, '-r');hold off
    legend(["data","model"])
end

function plot_attenuation(v, A, x, data)
    %A = v(1); % initial amplitude
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);

    figure;
    scatter(x, data, 4,'red'); hold on;
    plot(x, A_decay, '-b'); hold off
    legend(["data","model"])

end