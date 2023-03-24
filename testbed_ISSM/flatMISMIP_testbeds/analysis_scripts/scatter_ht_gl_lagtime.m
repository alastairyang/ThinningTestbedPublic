%% H(t) timeseries, terminus, and grounding line
% This study calculates the lag time between H(t) and GL/terminus timeseries 
% via cross-correlation

%% Main script
% parameters
md_type = 'expt'; % options: "expt", "ctrl"
%geom_type = 'deep'; % options: "deep", "shallow"
ds = 500; % meshgrid spacing

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
gl_lags = cell(2,n_simu);
front_lags = cell(2,n_simu);
dist_to_front = cell(2,n_simu);

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
        % correlation and cross-correlation
        % cross-correlation finds the delay time
    %     gl_corr = zeros(md.mesh.numberofvertices, 1);
        gl_lag = zeros(md.mesh.numberofvertices, 1);
    %     front_corr = zeros(md.mesh.numberofvertices, 1);
        front_lag = zeros(md.mesh.numberofvertices, 1);        
        dt = 0.1;
        % frontal retreat rate and grounding line retreat rate; smooth with
        % span == 10
        gl_rate = transpose(smooth(diff(ht_data.gl,1)/dt,10));
        front_rate = transpose(smooth(diff(ht_data.front,1)/dt,10));
        for k = 1:md.mesh.numberofvertices
            dhdt = smooth(diff(h_mat(:,k),1)/dt,10);
    %         % correlation
    %         coefs = corrcoef(gl_rate, dhdt);
    %         gl_corr(k) = coefs(1,2);
    %         coefs = corrcoef(front_rate, dhdt);
    %         front_corr(k) = coefs(1,2);
            % cross-correlation: we find the average of lags where correlations
            % is > 0.9
            % first, gl
            [c,lags] = xcorr(dhdt,gl_rate,'coeff');
            gl_lag(k) = nanmean(lags(c>0.9));
            % then, front
            [c,lags] = xcorr(dhdt,front_rate,'coeff');
            front_lag(k) = nanmean(lags(c>0.9));
        end
    
        [gl_lag_grid, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, gl_lag, ds);
        [front_lag_grid, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, front_lag, ds);
        % get the grounding line position at the last timestep
        ocean_mask = md.results.TransientSolution(end).MaskOceanLevelset;
        gl_dist = locate_groundingline(md, ocean_mask); 
        % extract center line, keep data behind the last grounding line
        yi_mid = floor(size(gl_lag_grid,1)/2);
        gl_lag_mid = gl_lag_grid(yi_mid, x<=gl_dist);
        front_lag_mid = front_lag_grid(yi_mid, x<=gl_dist);
    
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
        dist_to_front{q,j} = abs(front_x_mid - x_keep);

        % save to cells
        gl_lags{q,j} = gl_lag_mid;
        front_lags{q,j} = front_lag_mid;
    
        disp([modelname,' is processed.'])
    end
end

%% Making the plot
figure('Position',[100,100,800,800])
tiledlayout(3,3,'TileSpacing','none')
for j = 1:n_simu
    nexttile % plot in km and year
    % plot shallow
    scatter(dist_to_front{1,j}/1000, front_lags{1,j}/10,...
        'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on;
    scatter(dist_to_front{1,j}/1000, gl_lags{1,j}/10,...
        'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
    % plot deep
    scatter(dist_to_front{2,j}/1000, front_lags{2,j}/10,...
        '+','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on;
    scatter(dist_to_front{2,j}/1000, gl_lags{2,j}/10,...
        '+','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
    ylim([0,12])
    set(gca,'ytick',[1,5,9,13])
    set(gca,'xtick',[10,30,50])
%     if j == 1 % add legend
%         legend
%     end
%     if ismember(j, [1,4,7])
%         set(gca,'ytick',[0.5,0.7,0.9]);
%     else
%         set(gca,'ytick',[]);
%     end
%     if ismember(j, [7,8,9])
%         set(gca,'xtick',[1e4,3e4,5e4]);
%     else
%         set(gca,'xtick',[]);
%     end
end
ax = nexttile(4); ax.YLabel.String = 'Lag year';
ax = nexttile(8); ax.XLabel.String = 'Distance to calving front (m)';

exportgraphics(gcf,['plots/correlation_',md_type,'.pdf'],"ContentType","vector")
