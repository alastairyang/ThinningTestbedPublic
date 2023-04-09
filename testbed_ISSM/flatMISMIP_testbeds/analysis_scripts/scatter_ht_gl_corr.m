%% H(t) timeseries, terminus, and grounding line
% This study calculates the correlation between H(t) and GL/terminus timeseries 

%% Main script
% parameters
%md_type = 'expt'; % options: "expt", "ctrl"
ds = 100; % meshgrid spacing

% input
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
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

% pre-allocate
n_simu = size(folder_dir_groups{1}, 1);
gl_corrs = cell(2,2,n_simu);
front_corrs = cell(2,2,n_simu);
dist_to_front = cell(2,2,n_simu);
% gl_corrcoef_mean = zeros(2,2,n_simu);
% front_corrcoef_mean = zeros(2,2,n_simu);
%gl_decay_length = zeros(2,n_simu);
%front_decay_length = zeros(2,n_simu);

md_types = ["ctrl","expt"]; % want to plot both control and experiment in one plot
for tp = 1:length(md_types)
    md_type = md_types(tp);
    for q = 1:length(GLs) % shallow, deep
        % start extracting data
        for j = 1:n_simu
            % read the model
            group = folder_dir_groups{q};
            % load both the model data and extracted centerline data
            switch md_type
                case 'ctrl'
                    load([group.folder{j},'/', group.name{j}, '/', ctrl_name])
                    modelname = md.miscellaneous.name;
                    load(['analyzed_data/calve_only/ht_calve_model_',modelname(9:end),'.mat'])
                case 'expt' % here we use effective pressure experiment by default
                    load([group.folder{j},'/', group.name{j}, '/', expt_name])
                    modelname = md.miscellaneous.name;
                    load(['analyzed_data/mu_calve/ht_mu_calve_model_',modelname(9:end),'.mat'])
                otherwise 
                    warning('unsupported input')
            end
            [W, GL, FC] = parse_modelname(modelname);
            
            % put thickness change overtime into one big matrix
            result_tbl = struct2table(md.results.TransientSolution);
            n_time = size(result_tbl,1);
            h_mat = [];
            for i = 1:n_time
                h_mat = [h_mat; transpose(cell2mat(result_tbl.Thickness(i)))];
            end
    
            gl_corr = zeros(md.mesh.numberofvertices, 1);
            front_corr = zeros(md.mesh.numberofvertices, 1);
            dt = 0.1;
            % frontal retreat rate and grounding line retreat rate
            gl_rate = transpose(smooth(diff(ht_data.gl,1)/dt,30));
            front_rate = transpose(smooth(diff(ht_data.front,1)/dt,30));
            for k = 1:md.mesh.numberofvertices
                %dhdt = smooth(diff(h_mat(:,k),1)/dt,10);
                ht = h_mat(:,k);
                % correlation
                %coefs = corrcoef(gl_rate, dhdt);
                coefs = corrcoef(ht_data.gl, ht);
                gl_corr(k) = coefs(1,2);
                %coefs = corrcoef(front_rate, dhdt);
                coefs = corrcoef(ht_data.front, ht);
                front_corr(k) = coefs(1,2);
            end
        
            % interpolate onto a regular meshgrid
            [gl_corr_grid, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, gl_corr, ds);
            [front_corr_grid, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, front_corr, ds);
            % get the grounding line position at the last timestep
            gl_mask = md.results.TransientSolution(end).MaskOceanLevelset;
            gl_dist = locate_calvingfront(md, gl_mask); 
            % extract center line, keep data behind the last calving front
            yi_mid = floor(size(gl_corr_grid,1)/2);
            gl_corr_mid = gl_corr_grid(yi_mid, x<=gl_dist);
            front_corr_mid = front_corr_grid(yi_mid, x<=gl_dist);
        
            % interpolate to find the x of front at the center flow line
            ice_mask = md.results.TransientSolution(end).MaskIceLevelset;
            front_xy = isoline(md, ice_mask,'value',0);
            % to do interpolation, we always crop out the ridges to avoid repeating
            % values
            wid_factor = 0.3;
            front_y_crop = front_xy.y < max(front_xy.y)/2+W*wid_factor &...
            front_xy.y > max(front_xy.y)/2-W*wid_factor;
            front_y = front_xy.y(front_y_crop);
            front_x = front_xy.x(front_y_crop);
            front_x_mid = interp1(front_y, front_x, max(y)/2);
            x_keep = x(x <= gl_dist);
            dist_to_front{tp,q,j} = abs(front_x_mid - x_keep);
    
            % save to cells
            gl_corrs{tp,q,j} = gl_corr_mid;
            front_corrs{tp,q,j} = front_corr_mid;
    
    %         % estimate the decay length scale 
    %         options = optimset('Display','notify','PlotFcns',@optimplotfval);
    %         decay_v = [5]; % decay length scale initial guess
    %         % we only use the positive data points for decay length scale
    %         % estimation
    %         gl_corr_mid_c = gl_corr_mid(gl_corr_mid > 0);
    %         front_corr_mid_c = front_corr_mid(front_corr_mid > 0);
    %         % new x axis
    %         gl_xx = dist_to_front{q,j}(gl_corr_mid > 0)/1000;
    %         front_xx = dist_to_front{q,j}(front_corr_mid > 0)/1000;
    %         gl_decay_length(q,j)    = fminsearch(@minimize_attenuation, [decay_v,max(gl_corr_mid)], options, gl_xx, gl_corr_mid_c);
    %         front_decay_length(q,j) = fminsearch(@minimize_attenuation, [decay_v,max(front_corr_mid)], options, front_xx, front_corr_mid_c);
    %         % make a plot
    %         plot_attenuation(decay_length, max(gl_corr_mid), dist_to_perturb, peaks_downstream)
    %         % calculate mean correlation coef (10 km near the front)
    %         front_idx = find(dist_to_front{q,j}<10000);
    %         gl_corrcoef_mean(q,j) = mean(gl_corrs{q,j}(front_idx));
    %         front_corrcoef_mean(q,j) = mean(front_corrs{q,j}(front_idx));
        
            disp(['Model ', modelname(9:end),' is processed.'])
        end
    end
end

% %% Attentuation length scale (exponential decay model fit)
% [decay_length, A] = fminsearch(@minimize_attenuation, [decay_v,max(gl_corr_mid)], options, gl_xx, gl_corr_mid_c);
% plot_attenuation([decay_length, A], fliplr(gl_xx), fliplr(gl_corr_mid_c))
%% Making the plot
figure('Position',[100,100,500,500])
tiledlayout(3,3,'TileSpacing','none')
line_styles = [":","-"];
for p = 1:length(md_types)
    for j = 1:n_simu
        nexttile(j) % plot in km and year
        % plot deep
        dist_in_km = dist_to_front{p,2,j}/1000;
        plot_xi = find(dist_in_km<40);
        plot_x = dist_in_km(plot_xi);
        plot(plot_x, front_corrs{p,2,j}(plot_xi),...
            'LineStyle',line_styles(p), 'color',[253,187,132]/255,'LineWidth',2); hold on;
        plot(plot_x, gl_corrs{p,2,j}(plot_xi),...
            'LineStyle',line_styles(p),'color',[166,189,219]/255,'LineWidth',2); hold on
        % add a dot to the grounding line
        gl_x = plot_x(end);
        front_val = front_corrs{p,2,j}(end);
        gl_val = gl_corrs{p,2,j}(end);
        scatter(gl_x, front_val, 30, [253,187,132]/255,'filled'); hold on;
        scatter(gl_x, gl_val, 30, [166,189,219]/255,'filled'); hold on;
        ylim([0.6,1])
        xlim([0,40])
        %legend(["front","GL"],'Location','southwest')
        set(gca,'ytick',[1,5,9,13])
        set(gca,'xtick',[20,40])
        if ismember(j, [1,4,7])
            set(gca,'ytick',[0.8, 1]);
        else
            set(gca,'ytick',[]);
        end
        if ismember(j, [7,8,9])
            set(gca,'xtick',[20,40]);
        else
            set(gca,'xtick',[]);
        end
    end
end
ax = nexttile(4); ax.YLabel.String = 'Correlation'; ax.YLabel.FontSize = 15;
ax = nexttile(8); ax.XLabel.String = 'Distance to calving front (km)'; ax.XLabel.FontSize = 15;

exportgraphics(gcf,'plots/correlation_deep.png',"Resolution",600)

% for the shallow glaciers: but they don't have grounding lines so kind of
% pointless
% figure('Position',[600,100,500,500])
% tiledlayout(3,3,'TileSpacing','none')
% for j = 1:n_simu
%     nexttile % plot in km and year
%     % plot shallow
%     scatter(dist_to_front{1,j}/1000, front_corrs{1,j},...
%         'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on;
%     scatter(dist_to_front{1,j}/1000, gl_corrs{1,j},...
%         'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
%     ylim([0.5,1])
%     xlim([0,50])
%     legend(["front","GL"],'Location','southwest')
%     set(gca,'ytick',[1,5,9,13])
%     set(gca,'xtick',[10,30,50])
% end
% ax = nexttile(4); ax.YLabel.String = 'Correlation';
% ax = nexttile(8); ax.XLabel.String = 'Distance to calving front (km)';
% 
% exportgraphics(gcf,['plots/correlation_',md_type,'_shallow.png'],"Resolution",500)


%% APPENDIEX: functions
function err = minimize_attenuation(v, x, data)
    beta = v(1); % decay length scale
    A = v(2); % magnitude
    A_decay = A.*exp(-beta.*x);
    err = sqrt(mean((A_decay - data).^2));
end

function plot_attenuation(v, x, data)
    beta = v(1); % decay length scale
    A = v(2); % magnitude
    A_decay = A.*exp(-beta.*x);

    figure;
    scatter(x, data, 4,'red'); hold on;
    plot(x, A_decay, '-b'); hold off
    legend(["data","model"])

end
