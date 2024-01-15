function [px, py, pmag] = drivingstress_from_results(md, time_i)
%DRIVINGSTRESS_FROM_RESULTS calculate the driving stress, but not from
%md.geometry; from md.results.TransientSolution
%
%   Input:
%       md: ISSM model class
%       time_i: time index / No. row in solution table.

    % taken from function "slope.m"
    if dimension(md.mesh)==2
	    numberofelements=md.mesh.numberofelements;
	    numberofnodes=md.mesh.numberofvertices;
	    index=md.mesh.elements;
	    x=md.mesh.x; y=md.mesh.y;
    else
	    numberofelements=md.mesh.numberofelements2d;
	    numberofnodes=md.mesh.numberofvertices2d;
	    index=md.mesh.elements2d;
	    x=md.mesh.x2d; y=md.mesh.y2d;
    end

	surf=md.results.TransientSolution(time_i).Surface;
    %compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
    [alpha, beta]=GetNodalFunctionsCoeff(index,x,y);
    
    summation=[1;1;1];
    sx=(surf(index).*alpha)*summation;
    sy=(surf(index).*beta)*summation;
    s=sqrt(sx.^2+sy.^2);

    % end of modified "slope.m"
    % now start modified "drivingstress.m"

    % Average thickness over elements
    thickness = md.results.TransientSolution(time_i).Thickness;
    
    thickness_bar=(thickness(md.mesh.elements(:,1))+thickness(md.mesh.elements(:,2))+thickness(md.mesh.elements(:,3)))/3;
    
    % get driving stress for x and y components
    px=-md.materials.rho_ice*md.constants.g*thickness_bar.*sx;
    py=-md.materials.rho_ice*md.constants.g*thickness_bar.*sy;
    pmag=sqrt(px.^2+py.^2);
    
    
end

