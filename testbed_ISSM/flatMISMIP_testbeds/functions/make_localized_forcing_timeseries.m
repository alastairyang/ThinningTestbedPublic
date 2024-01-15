function [pulse_full, diffu_full, gauss_t] = make_localized_forcing_timeseries()
%CONSTRUCT_LOCALIZED_FORCING Make transient pulse and/or diffused pulse.
%This file assumes current working directory at flatMISMIP_testbeds

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

end

