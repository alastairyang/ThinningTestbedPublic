%% get distance to ice front for each grid point
ice_mask = md.results.TransientSolution(end).MaskIceLevelset;
front_xy = isoline(md, ice_mask,'value',0);
front_y_crop = front_xy.y < max(front_xy.y)/2+W*wid_factor &...
    front_xy.y > max(front_xy.y)/2-W*wid_factor;
front_y = front_xy.y(front_y_crop);
front_x = front_xy.x(front_y_crop);
front_x_interp = interp1(front_y, front_x, transpose(y(y_crop)));
front_x_interp = fillmissing(front_x_interp,"nearest");
front_xy_interp = [front_x_interp, transpose(y(y_crop))];
% distance
[X_crop, Y_crop] = meshgrid(x(x_crop),y(y_crop));
dist_to_front = abs(front_x_interp - X_crop);
