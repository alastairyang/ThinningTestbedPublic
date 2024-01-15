function C_update = mass_unloading(md, Hi, H0, k, C0, C_last, ocean_mask, m)
%MASS_UNLOADING calculate the change of coefficients in Weertman's sliding law
% as a result of change in ice thickness in the previous timestep
    
    g = md.constants.g;
    rho_ice = md.materials.rho_ice;
    
    % calculate new C-squared
    % if negative, real() -> 0 
    C_update = real(sqrt(C0.^2 + k.^2.*((rho_ice*g.*Hi).^(1/m)-...
                                        (rho_ice*g.*H0).^(1/m) ...
                                       ) ...
                         ) ...
                    );
    C_update(C_update > C0) = C0(C_update > C0); % avoid blowing up
end

