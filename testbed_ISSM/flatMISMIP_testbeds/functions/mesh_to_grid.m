function [grid_data, x, y] = mesh_to_grid(index, mesh_x, mesh_y, mesh_data, ds)
%MESH_TO_GRID convert meshed data to gridded data
%
%   Input:
%       index: md.mesh.elements
%       mesh_x: md.mesh.x
%       mesh_y: md.mesh.y
%       mesh_data: data of interest
%       ds: spatial interval
%
%   Output:
%       grid_data: gridded data
%       x: 1d vector of x axis
%       y: 1d vector of y axis


        Lx = max(mesh_x);
        Ly = max(mesh_y);
        x = 0:ds:Lx-ds;
        y = 0:ds:Ly-ds;
        grid_data = InterpFromMeshToGrid(index, mesh_x, mesh_y,...
                                         mesh_data, x, y, NaN);
end

