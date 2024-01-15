 function line_order = colormap_to_colororder(map, n, uplim, botskip)
%COLORMAP_TO_COLORORDER This function samples a colormap array (256x3 rgb)
%to multiple rgb units for line plot through interpolation.
%
%   input:
%       map [double array]: Nx3 rgb array (N is usually 256)
%       n [int]: number of lines we want to assign unique rgb to
%       uplim [int]: number of lines at the top of the map array we want to skip
%       botlim[int]: humber of lines at the bottom of the map array we want
%                    to skip
%
%   Output:
%       line_order [double array]: an array of size nx3

    if nargin<=2
        uplim = 1; %start from the beginning
        botskip = 0; % no skip
    end
    botlim = size(map,1)-botskip;
    row_axis = uplim:botlim; 
    row_interp_axis = linspace(uplim, botlim, n);
    row_interp_axis = floor(row_interp_axis);
    map_c = map(row_axis,:);
    
    line_order = interp1(row_axis, map_c, row_interp_axis);
    
end

