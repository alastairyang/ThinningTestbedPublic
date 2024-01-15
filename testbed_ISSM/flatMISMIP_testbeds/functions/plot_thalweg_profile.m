function [surface_profile, base_profile, x, y] = plot_thalweg_profile(md, ds, idx, make_plot)
%PLOT_PROFILE take the field data, interp to grid, and plot the first profile of
%values along the thalweg
%
%   Input:
%       md: ISSM model class
%       ds: spacing on a regular grid (in meter)
%       idx: numeral indices of solution in md.results.TransientSolution
%       make_plot: boolean, make a plot or not
%
%   Output:
%       data: 1d vector of requested variable along the thalweg
    
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
    x = 0:ds:Lx-ds;
    y = 0:ds:Ly-ds;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end

    % plot bed
    bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);

    % plot surface and base
    surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
        md.results.TransientSolution(idx).Surface,...
        x, y, NaN);
    base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
        md.results.TransientSolution(idx).Base,...
        x, y, NaN);

%     if nargin > 4 && ~isempty(mask)
%         switch mask
%             case 'ice'
%                 mask_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
%                     md.results.TransientSolution(mask_idx).MaskIceLevelset,...
%                     x, y, NaN);
%                 surface(mask_grid>0) = nan; % make ocean NaN
%                 base(mask_grid>0) = nan;
%             case 'grounded'
%                 mask_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
%                     md.results.TransientSolution(mask_idx).MaskOceanLevelset,...
%                     x, y, NaN);
%                 surface(mask_grid<0) = nan; % make floating ice NaN
%                 base(mask_grid<0) = nan;
%             otherwise
%                 fprintf('Unknown mask!\n')
%         end
%     end

    % profiles for surface and base for plotting
    surface_profile  = surface(mid_i,:);
    base_profile = base(mid_i,:);
    
    if nargin > 3 && make_plot
    figure('Position',[100,100,1000,300])
        plot(x, bed_profile, 'k');hold on;
        plot(x, surface_profile, 'Color','r'); hold on
        plot(x, base_profile,'Color','b');hold on
        ylim([-800,1800])
        time = md.results.TransientSolution(idx).time;
        title([md.miscellaneous.name, ', time = ', num2str(time)])
        hold off
    end
end

