function dH = plot_select_dhdt(md, xs, ys, make_plot)
%PLOT_SELECT_DHDT plot the h(t) or dh/dt(t) at selected ground control
%points

    nt = size(md.results.TransientSolution,2);
    % save the mesh elements and (x,y)
    mesh_x = md.mesh.x;
    mesh_y = md.mesh.y;
    % find grid index where we want to sample h(t)
    [ki, ~] = dsearchn([mesh_x, mesh_y], [xs, ys]);

    % initialize space to stare h(t)
    dH = zeros(length(ki), nt);

    % get all h(t) data
    for k = 1:length(ki)
        for j = 1:nt
            H = md.results.TransientSolution(j).Thickness;
            dH(k,j) = H(ki);
        end
    end
    % dH: change in H from the beginning
    dH = dH - dH(:,1);

    % time vector
    results = struct2table(md.results.TransientSolution);
    time = results.time;
    time = time - time(1);

    % figure
    if nargin > 3 && make_plot 
        color_length = n_point;
        red = [255, 51, 153]/255;
        sth = [153, 153, 255]/255;
        colors_p = [linspace(red(1),sth(1),color_length)',...
            linspace(red(2),sth(2),color_length)',...
            linspace(red(3),sth(3),color_length)'];
        figure;
        plot(time, dH);
        colororder(colors_p);
    end
    
end

