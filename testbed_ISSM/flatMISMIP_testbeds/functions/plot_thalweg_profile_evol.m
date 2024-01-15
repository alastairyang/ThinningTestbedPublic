function plot_thalweg_profile_evol(md,skip)
%PLOT_PROFILE take the field data, interp to grid, and plot the profile evolution of
%values along the thalweg

    if nargin == 1; skip = 1; end
    
    nt = size(md.results.TransientSolution,2);
    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)

    % color gradient
    color_length = nt;
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
                linspace(red(2),sth(2),color_length)',...
                linspace(red(3),sth(3),color_length)'];
    
    % geometry parameters
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 250; % spacing, 250 meter
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);

    figure
    bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);
    plot(x, bed_profile, 'black');hold on;
    for i = 1:skip:nt
        
        surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Surface,...
            x, y, NaN);
        base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Base,...
            x, y, NaN);
        bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
        surface_profile  = surface(mid_i,:);
        base_profile = base(mid_i,:);
        
        plot(x, bed_profile, 'k');hold on;
        plot(x, surface_profile, 'Color',colors_p(nt+1-i,:)); hold on
        plot(x, base_profile,'Color',colors_p(nt+1-i,:));hold on
        
        %hold off
        ylim([-800,1800])
        time = md.results.TransientSolution(i).time;
        title([md.miscellaneous.name, ', time = ', num2str(time)])
    end
    hold off
end

