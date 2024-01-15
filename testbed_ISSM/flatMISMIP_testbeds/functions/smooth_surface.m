function [new_data, X, Y] = smooth_surface(surface_data, new_x, new_y, n)
%SMOOTH_SURFACE This function smooths a surface
%
%   example: smooth_surface(surface_data, [1:10], [35:100])
%   
%   Input:
%       n: number of paddings added to each dimension
%   Output:
    
    [X, Y] = meshgrid(new_x, new_y);
    
    % add padding to the array as surface smoothing might bring noise to
    % the boundaries
    pad_data = padarray(surface_data, [n,n], 'replicate','both');
    % our smoothing patch is the same size as our padding
    smooth_data = medfilt2(pad_data, [n, n]);
    % remove the padding and restore the original dimension
    new_data = smooth_data(n+1:end-n, n+1:end-n);
    
    % visualize the data
    mesh(X, Y, new_data)
    
end

