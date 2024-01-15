function [dhdt, t] = sample_dhdt(md, xs)
%SAMPLE_DHDT Sample dh/dt timeseries at given point along the centerline
%
%   example: [dhdt, t] = sample_dhdt(md, [1e5,2e5,3e5])   
%
%   Input:
%       md: model class
%           ISSM model
%       xs: 1d array
%           distances from influx boundary (x=0) along the centerline.
%   
%   Output:
%       dhdt: 2d array
%           mxn array where m is number of points supplied in xs. m is the
%           number of simulation timestep. dh/dt
%       t: 1d array
%           simulation timesteps

    % build coordinate and find thalweg
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 250; % spacing, 250 meter
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        y_i = size(X,1)/2;
    else
        y_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(y_i,:);

    % get dh/dt at different points
    n_time = size(md.results.TransientSolution,2);
    n_points = length(xs);
    dhdt = zeros(n_points, n_time-1);

    for j = 1:length(xs)
        % find the nearest coordinate to the given x
        x_i = interp1(thalweg_x,1:length(thalweg_x), xs(j),'nearest');
        for i = 1:n_time
            if i == 1
                continue
            else
                dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
                surface_now  = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y, md.results.TransientSolution(i).Surface,x, y, NaN);
                surface_last = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,md.results.TransientSolution(i-1).Surface,x, y, NaN);
                dh = surface_now - surface_last;
                s_now  = surface_now(y_i, x_i);
                s_last = surface_last(y_i, x_i);
                
                dhdt(j,i) = (s_now - s_last)/dt;
            end
        end
    end
    results_tbl = struct2table(md.results.TransientSolution);
    t = results_tbl.time(1:end);
    
end

