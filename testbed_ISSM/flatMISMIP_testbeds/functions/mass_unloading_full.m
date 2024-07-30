function C_update = mass_unloading_full(md, Hi, H0, k, C0, pw, m)
%MASS_UNLOADING_FULL calculate the change of coefficients in Weertman's sliding law
% as a result of change in ice thickness in the previous timestep. We
% consider the full equation in the paper (B.14) which does not assume zero
% subglacial water pressure. Here the water pressure is assumed positive.
    
    g = md.constants.g;
    rho_ice = md.materials.rho_ice;
    
    % calculate new C-squared
    % if negative, real() -> 0 
    C_update = real(sqrt(C0.^2 + k.^2.*((rho_ice*g.*Hi - pw).^(1/m)-...
                                        (rho_ice*g.*H0 - pw).^(1/m) ...
                                       ) ...
                         ) ...
                    );
    C_update(C_update > C0) = C0(C_update > C0); % avoid blowing up
end

