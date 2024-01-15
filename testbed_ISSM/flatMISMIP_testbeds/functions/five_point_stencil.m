function data_d_full = five_point_stencil(data, h, dim)
%FIVE_POINT_STENCIL Finite difference stencil that uses five points, O(h^4)
%accuracy

    if dim == 1
        data = data';
    end
    % shift
    data_left1 = [data(:,end), data(:,1:end-1)]; % add 1 padding col to left (i+1)
    data_left2 = [data(:,end-1:end),data(:,1:end-2)]; % add 2 padding col to left (i+2)
    data_right1 = [data(:,2:end), data(:,1)]; % add 1 padding col to right (i-1)
    data_right2 = [data(:,3:end), data(:,1:2)]; % add 2 padding col to right (i-2)

    % derivative
    data_d = (-data_right2 + 8*data_right1 - 8*data_left1 + data_left2)/(12*h);
    data_d_crop = data_d(:,3:end-2);

    % add 4 nan columns to the boundary
    % you can use higher order stencils if those columns matter to you
    data_d_full = [nan*ones(size(data,1),2), data_d_crop, nan*ones(size(data,1),2)];
    
    % transpose back
    if dim == 1
        data_d_full = data_d_full';
    end

end

