function x_dist = locate_calvingfront(md,mask)
%LOCATE_GROUNDINGLINE find the calving front distance to x = 0 from
%MaskIceLevelset
    
    % use isoline to find 0
    % find centerline and output the x value there
    contours = isoline(md,mask,'value',0);
    % notice that there can be multiple closed contours
    % now search for each x arary, if mid-width y-coord is present
    found_flag = 0;
    y_tol = 200; % toleratance 200 m

    % There can be multiple contour groups. Here we iterate over each and
    % find the first one that is less than y_tol away from the centerline
    % (now this can be faulty, but hopefully it captures the actual contour
    % that is the grounding line)
    contour_group = 1;
    while ~found_flag
        width_mid = (max(md.mesh.y)-min(md.mesh.y))/2;
        [width_err, width_mid_i] = min(abs(contours(contour_group).y - width_mid));
        if width_err > y_tol
            % go to the next contour group
            contour_group = contour_group + 1;
        else
            found_flag = 1;
        end
    end
    x_dist = contours(contour_group).x(width_mid_i);
end

