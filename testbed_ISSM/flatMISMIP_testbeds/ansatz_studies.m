%% use diffusion ansatz to fit the data 
% get the dh/dt and x data from simulation
n_gcp = 20;
[dhdt, pos] = plot_sampled_dhdt(md, n_gcp);
%% fit model to data
parms_all = zeros(n_gcp,4);
%parms_all = zeros(n_gcp,3);
for gcp_i = 1:n_gcp
    %gcp_i = 20;
    x = 60000 - pos(gcp_i);
    data = dhdt(:,gcp_i);
    beta_guess = max(data) - min(data);
    disp_guess = min(data);
    v = [2000, 1000000, disp_guess, beta_guess]; % initial guess
    options = optimset('PlotFcns',@optimplotfval);
    parms = fminsearch(@minimize_ansatz, v, options, x, data);
    
    % plot model fit and data
    plot_ansatz(parms, x, data);hold on;
    parms_all(gcp_i,:) = parms;
end

% Peclet number
dist_terminus = 60000 - pos;
Pe = parms(:,1).*dist_terminus'./parms_all(:,2);


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
