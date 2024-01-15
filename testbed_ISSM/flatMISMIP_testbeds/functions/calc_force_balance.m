function [driving_S_smt, basal_R_smt, longi_grad_smt, later_grad_smt, x, y, gl_x, front_x] = calc_force_balance(md, ti, smooth_L,ds)
%CALC_FORCE_BALANCE Estimate the force balance at a given time from the
%model results. The estimate data is in a meshgrid form. The force balance
%components are smoothed to account for the coupling length scale.
%
%   Input:
%       md      [model class]: ISSM model
%       ti      [int/double]: the number of iteration (~time)
%       smooth_L[double]: smooth length scale (meter)
%
%   Output:
%       driving_S [double array]: driving stress
%       basal_R   [double array]: basal resistive stress (shear stress)
%       longi_grad[double array]: longitudinal stress gradient
%       later_grad[double array]: lateral stress gradient

    % parameter
    n = 3; % rheology nonlinearity
    n_ds = smooth_L/ds;
    %compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
    index = md.mesh.elements;

    % basal stress
    if size(md.friction.C,2) == 1
        % no sliding law coefficient change
        bs = md.friction.C.^2.*md.results.TransientSolution(ti).Vel/md.constants.yts;
    else
        % mass unloading experiment
        bs = md.friction.C(1:end-1,ti).^2.*md.results.TransientSolution(ti).Vel/md.constants.yts;
    end
    bs_list = bs(index);
    basal_R = mean(bs_list,2);
    % lateral
    Rxx = md.materials.rheology_B.*md.results.TransientSolution(ti).StrainRateeffective.^(1/n-1).*(2*md.results.TransientSolution(ti).StrainRatexx +   md.results.TransientSolution(ti).StrainRateyy);
    Rxy = md.materials.rheology_B.*md.results.TransientSolution(ti).StrainRateeffective.^(1/n-1).*md.results.TransientSolution(ti).StrainRatexy;
    % thickness
    H = md.results.TransientSolution(ti).Thickness;

    % driving stress
    driving_S = drivingstress_from_results(md, ti);

    % interp onto grids
    [driving_S, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, driving_S, ds);
    [basal_R, ~, ~]   = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, basal_R, ds);
    % smooth (along-flow direction) the driving stress at length scale
    % of 1 km
    % first we interpolate to mesh, then use 5-point FD stencil to get the
    % derivative
    [Rxxgrid,~,~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, Rxx, ds);
    [Rxygrid,~,~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, Rxy, ds);
    [Hgrid,x,y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, H, ds);
    longi_grad = -1*five_point_stencil(Rxxgrid.*Hgrid, ds, 2);
    later_grad = -1*five_point_stencil(Rxygrid.*Hgrid, ds, 1);

    % smooth: 
    % 1st element applies to smoothing window across the flow; 
    % 2nd element applies to smoothing window along the flow
    driving_S_smt  = imgaussfilt(driving_S, [1,n_ds]);
    longi_grad_smt = imgaussfilt(longi_grad, [1,n_ds]);
    later_grad_smt = imgaussfilt(later_grad, [1,n_ds]);
    basal_R_smt    = imgaussfilt(basal_R, [1,n_ds]);
    % set basal resistive stress outside the ocean level set to 0
    [ocean_mask, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, md.results.TransientSolution(ti).MaskOceanLevelset, ds);
    basal_R_smt(ocean_mask<=0) = 0;

    % find the grounding line location at this timestep
    gl_x = locate_groundingline(md,md.results.TransientSolution(ti).MaskOceanLevelset);
    front_x = locate_groundingline(md,md.results.TransientSolution(ti).MaskIceLevelset);
    
end

