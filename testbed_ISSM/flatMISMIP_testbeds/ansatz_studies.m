%% use diffusion ansatz to fit the data 

x = 20000;
v = [1000, 10000, -20, 20]; % initial guess
data = dhdt(:,15);
options = optimset('PlotFcns',@optimplotfval);
parms = fminsearch(@minimize_ansatz, v, options, x, data);

% plot model fit and data
plot_ansatz(parms, x, data)
disp(parms)


%% Ansatz
function err = minimize_ansatz(v, x, data)
    u = v(1); % diffusive wave velocity
    D = v(2); % diffusivity
    b = v(3); % displacement
    beta = v(4); % multiplier

    t = 0:0.1:26; % 26 years
    H = b + beta*erf((x-u.*t)./(sqrt(4.*D.*t)));
    err = mean((data - H').^2);
end

function plot_ansatz(v, x, data)
    u = v(1); % diffusive wave velocity
    D = v(2); % diffusivity
    b = v(3); % displacement
    beta = v(4); % multiplier

    t = 0:0.1:26; % 26 years
    H = b + beta*erf((x-u.*t)./(sqrt(4.*D.*t)));
    figure;
    plot(t, data,'-.b');hold on
    plot(t, H, '-r');hold off
    legend(["data","model"])
end
