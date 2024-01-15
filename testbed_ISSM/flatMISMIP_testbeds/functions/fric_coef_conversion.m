function C_target = fric_coef_conversion(law_from, law_to, md, C, H, Zb, m)
%FRIC_COEF_CONVERSION convert friction coefficients between Budd's sliding
%law and Weertman's sliding law
    
    rho_ice = md.materials.rho_ice;
    rho_water = md.materials.rho_water;
    g = md.constants.g;

    if strcmp(law_from, 'Budd')||strcmp(law_to, 'Weertman')
        % bed if above sl, make zero in effective pressure
        Zb(Zb>0) = 0;
        N = rho_ice.*H + rho_water.*Zb;
        N = max(N, 0);
        C_target = real(sqrt(C.^2*(g.*N).^(1/m)));

    elseif strcmp(law_from, 'Weertman')||strcmp(law_to, 'Budd')
        % bed if above sl, make zero in effective pressure
        Zb(Zb>0) = 0.0;
        N = rho_ice*H + rho_water*Zb;
        N = max(N, 0);
        C_target = real(sqrt(C.^2./((g.*N).^(1/m))));
        % if N is zero (some floating section) -> inf
        %C_target(isinf(C_target)) = 0.0;
        C_target(~isfinite(C_target)) = 0.0;
    else
        error('Laws conversion not supported yet!')
    end

end

