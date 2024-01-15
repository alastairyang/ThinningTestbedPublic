function [dhdt, pos] = plot_sampled_dhdt(md,n_point)
%PLOT_SAMPLED_DHDT plot the h(t) or dh/dt(t) at a series of ground control
%points along the center flowline

    sample_interval = 1000; % meter; distance between two control points
    ds = 50; % grid spacing, meter
    ds_i = sample_interval/ds;

    nt = size(md.results.TransientSolution,2);
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    x = 0:ds:Lx-ds;
    y = 0:ds:Ly-ds;
    [X,Y] = meshgrid(x, y);
    % find centerline index
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    % save the mesh elements and (x,y)
    mesh_x = md.mesh.x;
    mesh_y = md.mesh.y;
    % retrieve the last position of the ice front; spclevelset ~= 0
    pos = find(md.levelset.spclevelset(1:end-1,end) < 1 & md.levelset.spclevelset(1:end-1,end) > -1);
    x_front = min(mesh_x(pos)); %#ok<FNDSB> 
    % find grid index where we want to sample h(t)
    [~, x_i_nearest] = min(abs(x - x_front));
    % sampled points: start at 1 interval behind the last ice front
    %
    front_i = x_i_nearest-ds_i;
    end_i   = front_i - n_point*ds_i;
    sample_i = front_i:-ds_i:(end_i+ds_i);

    % remove model class; data store in table instead to clear space
    results = struct2table(md.results.TransientSolution);
    empty_md = md;
    empty_md.results.TransientSolution = [];
    clear md

    % initialize space to stare h(t)
    thalweg_sample_ht = [];

    % get all h(t) data
    for j = 1:nt
        surface   = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
            results.Thickness{j},...
            x, y, NaN);
        % elevations
        surface_profile  = surface(mid_i,sample_i);
        thalweg_sample_ht = [thalweg_sample_ht; surface_profile];
    end
    thalweg_sample_ht = thalweg_sample_ht - thalweg_sample_ht(1,:);

    % time vector
    time = results.time;
    time = time - time(1);

    % figure
    color_length = n_point;
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
        linspace(red(2),sth(2),color_length)',...
        linspace(red(3),sth(3),color_length)'];
    figure;
    plot(time, thalweg_sample_ht);
    colororder(colors_p);

    % construct x vectors where h(t) are taken
    pos = x(sample_i);
    % save to output
    dhdt = thalweg_sample_ht;

    
end

