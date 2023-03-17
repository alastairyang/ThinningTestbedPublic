%% Map view: analyze localized basal perturbation
% here we plot the map view of elevation change data and the relative
% grounding line movement

%% Experiment with polynomial de-trending
gauss_xloc = 3.2e4; % location of center of gaussian perturbation in meter
pulse_type = 'Diffu'; % types: "Diffu","Pulse"
geom_type = 'deep'; % types: "deep", "shallow"
expt_type = 'mu'; % types: "no_mu", "mu" (without mass-unloading; with mass unloading)

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);
% plot parameter 
ylabel_i = [1,4,7];
xlabel_i = [7,8,9];
Ws_symb = [10,20,30];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

% split the folder_dir into two groups, separated by grounding line depth
folder_dir_groups = cell(1,2);
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(foldernames_tbl,1),1);
    for j = 1:size(foldernames_tbl.name)
        GL_bool(j) = compare_GLvalue(foldernames_tbl.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = foldernames_tbl(find(GL_bool),:); %#ok<FNDSB> 
end

% we divide the dicussions by the grounding line depth
[~, shallowGL_i] = min(GLs);
[~, deeperGL_i]  = max(GLs);

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

switch expt_type
    case 'no_mu'
        ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
        expt_name = ['MISMIP_yangTransient_Calving_',pulse_type,'GaussianPerturb_8.mat'];
    case 'mu'
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
    otherwise
        error('Unknown experiment type!')
end

n_simu = size(folder_dir_groups{geom_i}, 1);
% % pre-allocate
W_symbs = zeros(n_simu,1);
FC_symbs = zeros(n_simu,3);
STs_cl = zeros(n_simu, 1200); % 1200 is length
LTs_cl = zeros(n_simu, 1200); 
% grounding line positions
gl_ctrl = zeros(n_simu,1);
gl_expt = zeros(n_simu,1);

figure('Position',[100,100,1000,600])
tiledlayout(2,3,'TileSpacing','none')
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % get model information
    modelname = md_ctrl.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname);
    % isolate the delta H from localized basal perturbation
    expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Surface{:}]', results_tbl_ctrl.time,'linear','extrap'));
    deltaH = expt_H_interp - [results_tbl_ctrl.Surface{:}];
    deltaH_cell = num2cell(deltaH,1);
    [md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
    % permute
    md_grid = permute(md_grid,[2,3,1]);
    
    % apply the last mask to data at all timesteps
    mask = mask_grid(end,:,:);
    for i = 1:size(md_grid,3)
        md_temp = md_grid(:,:,i);
        md_temp(mask >=0) = 0;
        md_grid(:,:,i) = md_temp;
    end

    % obtain the last grounding line positions from both the control and
    % experiment
    gl_ctrl(j) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(end).MaskOceanLevelset);
    gl_expt(j) = locate_groundingline(md_expt, md_expt.results.TransientSolution(end).MaskOceanLevelset);
    % timeseries decomposition
    % we reshape x,y into a long vector and after decomposed return to a map
    xl = size(md_grid,1); yl = size(md_grid,2); nt = size(md_grid,3);

    % Decompose the timeseries along the center flowline
    % only look at the perturbation years
    start_yr = 5;
    end_yr = 21;
    dt = 0.1;
    mid_i = floor(size(md_grid,1)/2);
    crop_ti = start_yr/dt:(end_yr/dt-1);
    sel_xi  = 1:20:size(md_grid,2);
    % crop out the extra time; select the equally spaced control points
    md_grid_mid = squeeze(md_grid(mid_i,sel_xi,:));
    md_grid_mid = md_grid_mid(:,crop_ti);
    LTs = zeros(size(md_grid_mid,1), length(crop_ti));
    % decomposition!
    for i = 1:size(md_grid_mid,1)
        [LT, ST] = trenddecomp(md_grid_mid(i,:));
        LTs(i,:) = LT;
    end
    STs = md_grid_mid - LTs;
    % color the lines by their distance to the gaussian perturbation
    % and save the image
    t_axis = linspace(start_yr, end_yr, size(LTs,2));
    sel_x = x(sel_xi);
    dist = sel_x - gauss_xloc; % positive: downstream
    abs_dist = abs(dist);
    unique_dist = unique(abs_dist);
    % first plot downstream ones
    color_length = length(unique_dist);
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = cool(color_length);
    colororder(colors_p)
    % color: from nearest to the furthest
    color_axis = transpose(linspace(min(abs_dist), max(abs_dist),length(dist))); % 10^0 to 10^4

    snps_plot_title = [md_ctrl.miscellaneous.name(9:end) '_' pulse_type '_' geom_type '_' expt_type '.png'];
%     figure('Position',[100,100,1100,500]);
%     % iterate over each line / control point
%     for i = 1:length(dist)
%         % Long term
%         subplot(2,2,[1,2])
%         line_color = interp1(color_axis, colors_p, abs_dist(i));
%         if dist(i) >= 0 % downstream
%             line_style = '-';
%         else % upstream
%             line_style = '--';
%         end
%         plot(t_axis, md_grid_mid(i,:),"Color",line_color,'LineStyle',line_style);
%         hold on
% 
%         subplot(2,2,3)
%         line_color = interp1(color_axis, colors_p, abs_dist(i));
%         if dist(i) >= 0 % downstream
%             line_style = '-';
%         else % upstream
%             line_style = '--';
%         end
%         plot(t_axis, LTs(i,:),"Color",line_color,'LineStyle',line_style);
%         hold on
% 
%         % short term
%         subplot(2,2,4)
%         plot(t_axis, STs(i,:),"Color",line_color,'LineStyle',line_style);
%         hold on
%     end
%     subplot(2,2,3);title('Long term'); xlabel('time');ylabel('m')
%     subplot(2,2,4);title('Short term'); xlabel('time');ylabel('m')
%     plot_title = [md_ctrl.miscellaneous.name(9:end),'_',pulse_type,'_',geom_type,'_',expt_type,'.png'];
%     exportgraphics(gcf,['plots/pulse_mu_plots/',plot_title],'Resolution',300)

%     % plotting animation: H(t)
%     [H_grid, x, y] = mesh_to_grid_overtime(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, num2cell(expt_H_interp,1), 50);
%     figure('Position',[100,100,1100,450]);
%     for i = 1:260
%         imagesc(x,y,squeeze(md_grid(:,:,i)));title(num2str(i/10));hold on;clim([-5,5]);colormap(diverg_colormap(50));colorbar;
%         pause(0.005);
%     end
%     % plotting animation: dH(t)/dt
%     dt = 0.1;
%     figure('Position',[100,100,1100,450]);
%     for i = 2:260
%         imagesc(x,y,squeeze((md_grid(:,:,i)-md_grid(:,:,i-1))/dt));title(num2str(i/10));hold on;clim([-5,5]);colormap(diverg_colormap(50));colorbar;
%         pause(0.005);
%     end

    is = [178,188];
%     if j == 3
%         is = [118, 198];
%     end
%     if j == 9
%         is = [118, 192];
%     end

    % retrieve grounding line position over time
    t_expt = results_tbl_expt.time;
    t_ctrl = results_tbl_ctrl.time;
    gls_expt = zeros(size(t_expt));
    gls_ctrl = zeros(size(t_ctrl));
    for i = 1:260
        gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset);
        gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset);
    end
    % interpolate
    gls_expt_interp = interp1(t_expt, gls_expt, t_ctrl);
    gls_diff = gls_expt_interp - gls_ctrl;

    % make static snapshot plots
%     datas = md_grid(:,:,is);
%     snapshots_fig = plot_result_snapshots(x, y, datas, W);
% 
%     % PLOT!
%     nexttile(3,[2,1])
%     t_shift = t_ctrl-t_ctrl(1);
%     plot(gls_diff, t_shift,'-k','LineWidth',1); 
%     % add the times of the two snapshots
%     hold on;
%     scatter(gls_diff(is), t_shift(is),40,'filled','red');
%     set(gca, 'YDir','reverse'); ylim([0,26]);
%     xlabel('Relative grounding line position (m)','Interpreter','latex','FontSize',13)
%     ylabel('Time (yr)','Interpreter','latex','FontSize',13)
% 
%     % save the graphics
%     exportgraphics(gcf,['plots/pulse_mu_plots/snapshots/' snps_plot_title], 'Resolution', 300)
    
    % make a GIF
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
    fig_length = 900;
    n_snapshots = 1;
    fig_width  = 400;

    gif(['plots/gifs/',snps_plot_title(1:end-4),'.gif'])
    fig = figure('Position',[100,100,fig_length, fig_width]);
    tiledlayout(2,4,'TileSpacing','compact')
    for ii = 1:260
        % first plot the map view of H(t)
        gif('frame',gcf)
        nexttile(1,[1,2])
        data = squeeze(md_grid(:,:,ii));
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
        colorbar

        % then plot the map view of dH(t)/dt
        nexttile(5,[1,2])
        dt = 0.1;
        if ii == 1
            data = zeros(size(md_grid(:,:,ii)));
        else
            data = squeeze(md_grid(:,:,ii)-md_grid(:,:,ii-1))/dt;
        end
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
        colorbar

        % add grounding line
        nexttile(3,[2,1])
        t_shift = t_ctrl-t_ctrl(1);
        plot(gls_diff, t_shift,'-k','LineWidth',1); 
        % add the times of the two snapshots
        hold on;
        % add the time as a dot
        scatter(gls_diff(ii), t_shift(ii),40,'filled','red');
        set(gca, 'YDir','reverse'); ylim([0,26]);
        xlabel('Relative grounding line position (m)','Interpreter','latex','FontSize',13)
        ylabel('Time (yr)','Interpreter','latex','FontSize',13)
        pause(0.05)
        sgtitle(['Time = ',num2str(ii/10)]);

        % add forcing plot
        nexttile(4,[2,1])
        switch pulse_type
            case "Diffu"
                [~, pulse, gauss_t] = make_localized_forcing_timeseries();
            case "Pulse"
                [pulse, ~, gauss_t] = make_localized_forcing_timeseries();
            otherwise
                error('Unknown pulse type!')
        end
        plot(pulse, gauss_t,'-k','LineWidth',1); 
        % add the times of the two snapshots
        hold on;
        % add the time as a dot
        pulse_dot = interp1(gauss_t, pulse, t_shift(ii));
        scatter(pulse_dot, t_shift(ii),40,'filled','red'); hold off
        set(gca, 'YDir','reverse'); ylim([0,26]);
        xlabel('$\alpha$','Interpreter','latex','FontSize',13)
        ylabel('Time (yr)','Interpreter','latex','FontSize',13)

        pause(0.05)
        sgtitle(['Time = ',num2str(ii/10)]);

    end
    
    % report
    disp(['model ',snps_plot_title(1:end-4), ' is processed!'])
end

%% OUTDATED: Polynomial Detrending iteration
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % get model information
    modelname = md_ctrl.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname);
    % isolate the delta H from localized basal perturbation
    expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
    deltaH = expt_H_interp - [results_tbl_ctrl.Thickness{:}];
    deltaMask = [results_tbl_expt.MaskOceanLevelset{:}] - [results_tbl_ctrl.MaskOceanLevelset{:}];
    deltaH_cell = num2cell(deltaH,1);
    [md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
    % permute
    md_grid = permute(md_grid,[2,3,1]);
    
    % apply the last mask to data at all timesteps
    mask = mask_grid(end,:,:);
    for i = 1:size(md_grid,3)
        md_temp = md_grid(:,:,i);
        md_temp(mask >=0) = 0;
        md_grid(:,:,i) = md_temp;
    end

    % obtain the grounding line positions from both the control and
    % experiment
    gl_ctrl(j) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(end).MaskOceanLevelset);
    gl_expt(j) = locate_groundingline(md_expt, md_expt.results.TransientSolution(end).MaskOceanLevelset);
    % timeseries decomposition
    % we reshape x,y into a long vector and after decomposed return to a map
    xl = size(md_grid,1); yl = size(md_grid,2); nt = size(md_grid,3);
    md_grid_v = transpose(reshape(md_grid, [xl*yl, size(md_grid,3)]));


    STs = detrend(md_grid_v, 18); 
    LTs = md_grid_v - STs;
    LTs = reshape(LTs, [nt, xl, yl]);
    STs = reshape(STs, [nt, xl, yl]);
    % permute back (time axis at the 3rd dimension)
    LTs = permute(LTs, [2,3,1]);
    STs = permute(STs, [2,3,1]);
    
    % animate
%     figure('Position',[100,100,800,150]);for i = 1:nt; imagesc(x,y,squeeze(LTs(:,:,i)));clim([-10,0]);colorbar;pause(0.05);end 
%     figure;for i = 1:nt; imagesc(x,y,squeeze(STs(:,:,i)));clim([-10,0]);pause(0.05);end
    
    % find the spatial mean and spatial std
    LTs_crop = LTs(:,:,51:end); STs_crop = STs(:,:,51:end); % skip the first 5 years
    LTs_mean = mean(LTs_crop, 3); STs_max = max(STs_crop,[], 3);STs_min = min(STs_crop,[], 3);
    LTs_std  = std(LTs_crop, 0, 3);  STs_std  = std(STs_crop, 0, 3);
    % long term mean timeseries
    LTs_time_mean = squeeze(nanmean(LTs, [1,2]));
    LTs_time_min  = squeeze(min(LTs, [],[1,2]));
    
    md_name = md_ctrl.miscellaneous.name(9:end);
    % trend: find total delta H
    % cyclic: find range (max - min)
    last_LTs = -1*squeeze(LTs(:,:,end));
    STs_range = STs_max - STs_min;

    % get the corresponding symbols for this scatter plot
    W_symbs(j,1) = Ws_symb(W==Ws);
    FC_symbs(j,:) = FCs_symb(FC==FCs,:);
    % save the plot
%     save_dir = ['plots/diffu_mu_plots/',md_name,'.png'];
%     exportgraphics(gcf, save_dir, 'Resolution',300)
    disp(['model ',md_name,' is completed!'])

    % save the center flow line
    STs_cl(j,:) = STs_range(size(STs_range,1)/2,:);
    LTs_cl(j,:) = last_LTs(size(last_LTs,1)/2,:);

    % report
    disp(['model ', modelname(9:end),' type ',geom_type, ' is processed!'])
end
%% OUTDATED: accompanied plot for detrending method: plot the along flow profile
ds = 50;
ylim_up = 25;
% parameter for gaussian patch dimensions
gauss_xloc = 3.2e4;
gauss_width_ratio = 0.08;
wid_eff_factor = 4; % 4 times the sigma (in the gaussian width) 

figure('Position',[100,100,800,600])
tiledlayout(2,1,'TileSpacing','none')
% cyclic component
nexttile
for i = 1:n_simu
    % make a new along flow x axis that have the x positions of the
    % grounding line
    old_x = 0:ds:size(STs_cl,2)*ds-1;
    new_x = sort([old_x, gl_ctrl(i), gl_expt(i)]);
    STs_cl_interp = interp1(old_x, STs_cl(i,:), new_x);
    % plot solid line from x = 0 to gl from the experiment
    gl_expt_i = find(new_x==gl_expt(i));
    plot(new_x(1:gl_expt_i), STs_cl_interp(1:gl_expt_i),'LineStyle','-','LineWidth',W_symbs(i,1)*0.05,'Color',FC_symbs(i,:));
    hold on
    plot(new_x(gl_expt_i+1:end), STs_cl_interp(gl_expt_i+1:end),'LineStyle','-.','LineWidth',W_symbs(i,1)*0.05,'Color',FC_symbs(i,:));
    hold on
    % add the gl position as points
    scatter(new_x(new_x==gl_expt(i)), STs_cl_interp(new_x==gl_expt(i)),W_symbs(i,1)*1.5,...
            'filled','o','MarkerFaceColor',FC_symbs(i,:));hold on
    scatter(new_x(new_x==gl_ctrl(i)), STs_cl_interp(new_x==gl_ctrl(i)),W_symbs(i,1)*1.5,...
            'o','MarkerEdgeColor',FC_symbs(i,:));hold on
end
hold off
ylim([0,ylim_up])
for j = 1:length(Ws)
    patch_width = (Ws(j)*gauss_width_ratio*wid_eff_factor)/2;
    left_x = gauss_xloc - patch_width;
    right_x = gauss_xloc + patch_width;
    pt = patch([left_x,right_x,right_x,left_x],[0,0,ylim_up,ylim_up],[217,217,217]/255);
    pt.FaceAlpha = 0.3;
    pt.EdgeAlpha = 0.5;
end
text_x = 1;
text_y = ylim_up - 1;
text(text_x, text_y, 'Cyclic magnitude','Interpreter','latex','FontSize',13)
ylabel('Elevation change (m)','Interpreter','latex','FontSize',13)

% Trend component
nexttile
for i = 1:n_simu
    % make a new along flow x axis that have the x positions of the
    % grounding line
    old_x = 0:ds:size(LTs_cl,2)*ds-1;
    new_x = sort([old_x, gl_ctrl(i), gl_expt(i)]);
    LTs_cl_interp = interp1(old_x, LTs_cl(i,:), new_x);
    % plot solid line from x = 0 to gl from the experiment
    gl_expt_i = find(new_x==gl_expt(i));
    plot(new_x(1:gl_expt_i), LTs_cl_interp(1:gl_expt_i),'LineStyle','-','LineWidth',W_symbs(i,1)*0.05,'Color',FC_symbs(i,:));
    hold on
    plot(new_x(gl_expt_i+1:end), LTs_cl_interp(gl_expt_i+1:end),'LineStyle','-.','LineWidth',W_symbs(i,1)*0.05,'Color',FC_symbs(i,:));
    hold on
    % add the gl position as points
    scatter(new_x(new_x==gl_expt(i)), LTs_cl_interp(new_x==gl_expt(i)),W_symbs(i,1)*1.5,...
            'filled','o','MarkerFaceColor',FC_symbs(i,:));hold on
    scatter(new_x(new_x==gl_ctrl(i)), LTs_cl_interp(new_x==gl_ctrl(i)),W_symbs(i,1)*1.5,...
            'o','MarkerEdgeColor',FC_symbs(i,:));hold on
end
ylim([0,ylim_up])
hold off
for j = 1:length(Ws)
    patch_width = (Ws(j)*gauss_width_ratio*wid_eff_factor)/2;
    left_x = gauss_xloc - patch_width;
    right_x = gauss_xloc + patch_width;
    pt = patch([left_x,right_x,right_x,left_x],[0,0,ylim_up,ylim_up],[217,217,217]/255);
    pt.FaceAlpha = 0.3;
    pt.EdgeAlpha = 0.5;
end
xlabel('Distance (m)','Interpreter','latex','FontSize',13)
ylabel('Elevation change (m)','Interpreter','latex','FontSize',13)
text_x = 1;
text_y = ylim_up-1;
text(text_x, text_y, 'Total thinning in trend component','Interpreter','latex','FontSize',13);
plot_name = ['plots/local_basal_profile_', expt_type,'_',pulse_type, '_', geom_type, '.png'];
exportgraphics(gcf, plot_name ,'Resolution',300)


%% APPENDIX: un-used code
% codes below are explorations of other methods.
%% Experiment with EOF
% parameter
n_m = 5;

md_expt = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat').md;
md_ctrl = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
deltaH = expt_H_interp - [results_tbl_ctrl.Thickness{:}];
deltaMask = [results_tbl_expt.MaskOceanLevelset{:}] - [results_tbl_ctrl.MaskOceanLevelset{:}];
deltaH_cell = num2cell(deltaH,1);
[md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
[mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
% permute
md_grid = permute(md_grid,[2,3,1]);

% apply the last mask to data at all timesteps
mask = mask_grid(end,:,:);
for i = 1:size(md_grid,3)
    md_temp = md_grid(:,:,i);
    md_temp(mask >=0) = 0;
    md_grid(:,:,i) = md_temp;
end
    

% EOF
[eofmaps, pc, expvar] = eof(md_grid,5);

% Plot
figure; 
subplot(2,5,1); imagesc(x, y,eofmaps(:,:,1)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 1')
subplot(2,5,2); imagesc(x, y,eofmaps(:,:,2)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 2')
subplot(2,5,3); imagesc(x, y,eofmaps(:,:,3)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 3')
subplot(2,5,4); imagesc(x, y,eofmaps(:,:,4)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 4')
subplot(2,5,5); imagesc(x, y,eofmaps(:,:,5)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 5')
% temporal evolution
t = linspace(0,26,size(pc,2));
subplot(2,5,6); anomaly(t,pc(1,:));
subplot(2,5,7); anomaly(t,pc(2,:));
subplot(2,5,8); anomaly(t,pc(3,:));
subplot(2,5,9); anomaly(t,pc(4,:));
subplot(2,5,10); anomaly(t,pc(5,:));

%% Experiment with DMD
% use two models as an illustration (localized basal perturbation)
% parameter
n_m = 20;

md_expt = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat').md;
md_ctrl = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
expt_H_interp = transpose(interp1(results_tbl_expt.time, [results_tbl_expt.Thickness{:}]', results_tbl_ctrl.time,'linear','extrap'));
deltaH = expt_H_interp - [results_tbl_ctrl.Thickness{:}];
deltaMask = [results_tbl_expt.MaskOceanLevelset{:}] - [results_tbl_ctrl.MaskOceanLevelset{:}];
deltaH_cell = num2cell(deltaH,1);
[md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
[mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
% permute
md_grid = permute(md_grid,[2,3,1]);

% apply the last mask to data at all timesteps
mask = mask_grid(end,:,:);
for i = 1:size(md_grid,3)
    md_temp = md_grid(:,:,i);
    md_temp(mask >=0) = 0;
    md_grid(:,:,i) = md_temp;
end

% DMD
% make the first dimension the time dimension
md_grid = permute(md_grid, [3,1,2]);
% truncate in time
md_grid = md_grid(50:210,:,:);
[evals, evecs, mode_amps, mode_freq, growth_rates, mode_E] = dmd_rom(md_grid, n_m, 0.1); % first 5 mode
dominant_mode_No=find(mode_amps==max(mode_amps));
%==========Mode Amplitudes===========
figure('color', 'w');
stem(abs(mode_freq),abs(mode_amps), 'k^', 'filled'); hold on
plot(abs(mode_freq(dominant_mode_No)),abs(mode_amps(dominant_mode_No)),'ro','MarkerSize',10);
set(gca, 'YScale', 'log')
title('Mode amplitude versus Frequency');
xlabel('f [Hz]');ylabel('Amplitude');
%==========Evolve one mode in time============
m=19; % plot the first m mode
theta=linspace(0,16,160);
figure('color', 'w','Position',[100,100,1000,300]); 
for i=1:length(theta)
    H=squeeze(real(mode_amps(m)*evecs(m,:,:)*exp(1i*evals(m)*theta(i))));
    imagesc(x,y,H)
    title(['time is ',num2str(theta(i))])
    clim([-3, 3])
    colorbar
    pause(0.02);
end
% all first n modes
figure;
theta=linspace(0,6*pi,100);
for j = 1:5
    subplot(5,1,j)
    H=squeeze(real(evecs(j,:,:)*exp(1i*theta(3))));
    imagesc(x,y,H)
    drawnow;
end
% Now, how to reconstruct the sinusoidal wave?
theta=linspace(0,16,160);
m_sel = [13,14,15,16];
figure('color', 'w','Position',[100,100,1000,300]); 
for i=1:length(theta)
    H=squeeze(sum(real(mode_amps(m_sel).*evecs(m_sel,:,:).*exp(1i*evals(m_sel)*theta(i))),1));
    imagesc(x,y,H)
    title(['time is ',num2str(theta(i))])
    clim([-10, 10])
    colorbar
    pause(0.02);
end