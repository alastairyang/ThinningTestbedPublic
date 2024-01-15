function diffu_gauss = make_diffu_gauss(pulse_tscale, diffu_tscale, efold, total_tscale, dt)
%MAKE_DIFFU_GAUSS make the diffused gaussian pulse (temporal osccillation
%of the basal perturbation
    
    t = (-total_tscale/2):dt:(total_tscale/2 - dt);
    alpha = efold/((diffu_tscale/2)^2);
    diffu_gauss = (pulse_tscale/diffu_tscale)*exp(-alpha*t.^2);

end

