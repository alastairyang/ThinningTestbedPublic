function [grid_datas, x, y] = mesh_to_grid_overtime(index, mesh_x, mesh_y, data, ds)
%MESH_TO_GRID_OVERTIME This script interpolate field data from ISSM results
%over multiple timesteps into a grid data organized in a matrix. It is
%basically a multi-time-step version of 'mesh_to_grid.m' function. The
%input 'data' should be a 1-d cell array.

    Lx = max(mesh_x);
    Ly = max(mesh_y);
    x = 0:ds:Lx-ds;
    y = 0:ds:Ly-ds;
    n_t = length(data);
    
    grid_datas = zeros(n_t, length(y), length(x));
    for i = 1:n_t
        grid_datas(i,:,:) = InterpFromMeshToGrid(index, mesh_x, mesh_y,...
            data{i}, x, y, NaN);
    end


end

