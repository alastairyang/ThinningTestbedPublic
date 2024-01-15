function pulse_gauss = make_pulse_gauss(pulse_tscale, efold, total_tscale, dt, t_shift)
%MAKE_PULSE_GAUSS make a short gaussian pulse (temporal osccillation
%of the basal perturbation

    t = (-total_tscale/2):dt:(total_tscale/2 - dt);
    alpha = efold/((pulse_tscale/2)^2);
    pulse_gauss = exp(-alpha*(t+t_shift).^2);
end

