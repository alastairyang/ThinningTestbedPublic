function rgb = diverg_colormap(n)
%DIVERG_COLORMAP Source: https://blogs.mathworks.com/steve/2015/01/20/divergent-colormaps/

    rgb = [ ...
        94    79   162
        50   136   189
       102   194   165
       171   221   164
       230   245   152
       255   255   191
       254   224   139
       253   174    97
       244   109    67
       213    62    79
       158     1    66  ] / 255;
    
    cb_axis = 1:size(rgb,1);
    if n > size(rgb,1)
        new_axis = linspace(1, size(rgb,1), n);
        rgb = interp1(cb_axis, rgb, new_axis);
    end
    
    % illustrate the colormap
%     b = repmat(linspace(0,1,200),20,1);
%     imshow(b,[],'InitialMagnification','fit')
%     colormap(rgb)
end

