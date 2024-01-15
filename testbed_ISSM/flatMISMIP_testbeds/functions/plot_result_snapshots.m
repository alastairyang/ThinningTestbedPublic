function fig = plot_result_snapshots(x,y,datas,W)
%PLOT_RESULT_SNAPSHOTS
%
%   Input:
%       x
%       y
%       data [double array]: 3D data array; 3rd dimension is time (different snapshots)
%       W [double]: fjord width

    ds = mean(x(2:end)-x(1:end-1));
    W_eff = 0.5*W;
    x_up  = 1e4; % upstream limit for cropping
    x_down = 5e4; % downstream limit for cropping, close to the initial terminus
    % cropping parameters
    yi_mid = floor(length(y)/2);
    yi_low = yi_mid - floor(W_eff/ds);
    yi_up  = yi_mid + floor(W_eff/ds);
    xi_up = floor(x_up/ds);
    xi_down = floor(x_down/ds);

    n_snapshots = size(datas,3);
    fig_length = 900;
    fig_width  = n_snapshots*(150);
    fig = figure('Position',[100,100,fig_length, fig_width]);
    tiledlayout(n_snapshots,3,'TileSpacing','compact')
    for i = 1:n_snapshots
        %subplot(n_snapshots,1,i)
        nexttile(3*i-2,[1,2])
        data = squeeze(datas(:,:,i));
        % cropped axis
        data_c = data(yi_low:yi_up, xi_up:xi_down);
        x_c = x(xi_up:xi_down);
        y_c = y(yi_low:yi_up);
        
        % mesh plot
        s = imagesc(x_c, y_c, data_c); hold on;
        %view(30, 35)
        set(gca,'DataAspectRatio',[1 1 1]) % axis scaling; emphasize z axis
        axis off;
        grid off; box off; 
        colormap(diverg_colormap(50)); clim([-7,7]);
    end

end

