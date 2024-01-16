%% Plot the map view 2D velocity, thickness, basal drag at the steady state; 
% We also generate a table that summarizes the near-terminus averages of
% these observables.
ds = 50; % grid size for the regular grid
mid_y = 6000; % half width of the domain width (constant for all testbeds)
front_x = 56650; % terminus distance to x = 0

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
md_name = 'MISMIP_yangTransient_Steadystate_Extended.mat';
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);

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
vels = cell(2, n_simu);
hs = cell(2, n_simu);
tau_bs = cell(2, n_simu);
gls = cell(2, n_simu);
terminus = cell(2, n_simu);
float_L = cell(2, n_simu);
gl_depths = cell(2, n_simu);
width_eff = cell(2, n_simu);
% create a table to record model parameter
md_vartypes = ["int16","int16","double"];
md_varnames = ["W","GL","FC"];
md_size = [18,3];
md_params = table('Size',md_size, 'VariableTypes', md_vartypes, 'VariableNames', md_varnames);

% extrac data
for i = 1:length(GLs) % iterate over the two grounding line depths
    for j = 1:n_simu % iterate over different width and sliding law coef
        % read the model
        group = folder_dir_groups{i};
        md = load([group.folder{j},'/', group.name{j}, '/', md_name]).md;
        % query data from the last time step
        vel  = md.results.TransientSolution(end).Vel;
        h    = md.results.TransientSolution(end).Thickness;
        base = md.results.TransientSolution(end).Base;
        % get the sliding law coefficient and calculate the basal drag
        tau_b = md.friction.C.^2.*(vel/md.constants.yts);
        mask = md.results.TransientSolution(end).MaskOceanLevelset;
        tau_b(mask<0) = nan;
        
        % interp to regular grid
        [vels{i,j}, ~, ~]  = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, vel, ds);
        [hs{i,j}, ~, ~]    = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, h, ds);
        [tau_bs{i,j}, ~, ~]= mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, tau_b, ds);
        [base, x, y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, base, ds);
        [X, Y] = meshgrid(x,y);
        % get the grounding line and terminus positions and depths
        gls{i,j} = isoline(md, md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
        %terminus{i,j} = isoline(md, md.results.TransientSolution(end).MaskIceLevelset,'value',0);
        % if multiple isoline contours were found, we only save the one that has
        % most nodes (most likely the grounding line)
        if size(gls{i,j},2) ~= 1
            gl_tbl = struct2table(gls{i,j});
            gl_tbl = sortrows(gl_tbl, 'nods','descend');
            gls{i,j} = table2struct(gl_tbl(1,:));
        end

        modelname = md.miscellaneous.name(9:end);
        % save the model parameter
        l_idx = sub2ind([length(GLs),n_simu], i, j);
        [md_params.W(l_idx), md_params.GL(l_idx), md_params.FC(l_idx)] = parse_modelname(md.miscellaneous.name);
        % get the corresponding depth there
        % only look at the central flow trunk
        gl_yi_keep = find(gls{i,j}.y < mid_y + md_params.W(l_idx)/2*0.5 & gls{i,j}.y > mid_y - md_params.W(l_idx)/2*0.5);
        gl_y_keep = gls{i,j}.y(gl_yi_keep); 
        gl_x_keep = gls{i,j}.x(gl_yi_keep);
%         % similarly, for the terminus position
%         tm_yi_keep = find(terminus{i,j}.y < mid_y + md_params.W(l_idx)/2*0.5 & terminus{i,j}.y > mid_y - md_params.W(l_idx)/2*0.5);
%         tm_y_keep = terminus{i,j}.y(tm_yi_keep); 
%         tm_x_keep = terminus{i,j}.x(tm_yi_keep);
        
        % interp2 to get the grounding line at the sampled points
        [Xq, Yq] = meshgrid(gl_x_keep, gl_y_keep);
        gl_depths{i,j} = interp2(X, Y, base, gl_x_keep, gl_y_keep);

        % floating section length
        if gl_depths{i,j} > -200 % fully grounded
            float_L{i,j} = 0;
        else
            float_L{i,j} = front_x - gl_x_keep;
        end
        
        % effective width (0.7 of the prescribed length)
        width_eff{i,j} = md_params.W(l_idx)*0.7;
        
        disp(['Model ', modelname, ' processing is completed!'])
    end
end
% get the x,y axis
[~, x, y]= mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, tau_b, ds);
clear md

%% Create the table and export
% parameters
crop_length = 10000; % 10 km behind the terminus
% recording mean, max, and min values for each map variable
vartypes = ["double","double","double"];
varnames = ["vel","h","tau_b"]; % velocity, thickness, and basal drag
size = [18, length(vartypes)];
ss_mean_tbl = table('Size',size,'VariableTypes',vartypes,'VariableNames',varnames);
ss_min_tbl = table('Size',size,'VariableTypes',vartypes,'VariableNames',varnames);
ss_max_tbl = table('Size',size,'VariableTypes',vartypes,'VariableNames',varnames);
% record mean for other non-map observables 
vartypes = ["double","double","double"];
varnames = ["width_eff","gl_depth","float_length"]; % effective width, grounding line depth, and floating sectioin length
size = [18, length(vartypes)];
ss_tbl = table('Size',size,'VariableTypes',vartypes,'VariableNames',varnames);

for i = 1:length(GLs)
    for j = 1:n_simu
        % first get linear index
        l_idx = sub2ind([length(GLs),n_simu], i, j);
        W  = md_params.W(l_idx);
        % record map variables
        [ss_min_tbl.vel(l_idx),   ss_max_tbl.vel(l_idx),   ss_mean_tbl.vel(l_idx)]  = crop_domain_stats(vels{i,j}, W, crop_length, front_x, ds, x, y);
        [ss_min_tbl.h(l_idx),     ss_max_tbl.h(l_idx),     ss_mean_tbl.h(l_idx)]    = crop_domain_stats(hs{i,j}, W, crop_length, front_x, ds, x, y);
        [ss_min_tbl.tau_b(l_idx), ss_max_tbl.tau_b(l_idx), ss_mean_tbl.tau_b(l_idx)] = crop_domain_stats(tau_bs{i,j}, W, crop_length, front_x, ds, x, y);
        % record non-map variables 
        ss_tbl.width_eff(l_idx)    = mean(width_eff{i,j},'all');
        ss_tbl.gl_depth(l_idx)     = mean(gl_depths{i,j},'all');
        ss_tbl.float_length(l_idx) = mean(float_L{i,j},'all');
    end
end
% export to spreadsheets
writetable(ss_tbl, 'result_tables/ss_md_params.csv');
writetable(ss_min_tbl, 'result_tables/ss_md_field_mins.csv')
writetable(ss_max_tbl, 'result_tables/ss_md_field_maxs.csv')
writetable(ss_mean_tbl, 'result_tables/ss_md_field_means.csv')

%% make a tiled plot
%% velocity
fig_name = 'vel.svg';
f = figure('Position',[100,100,1000,500]);

p1=uipanel('parent',f);
p2=uipanel('parent',f);
p1.Title=("Shallow"); p1.BackgroundColor = 'w';
p2.Title=("Deep");    p2.BackgroundColor = 'w';
p1.Position=[0,0.5,1,0.5]; 
p2.Position=[0,0,1,0.5];
T1=tiledlayout(p1,3,3,'TileSpacing','none');
T2=tiledlayout(p2,3,3,'TileSpacing','none');
% colorbar range
cb_top = 4; cb_btm = 0;

plot_md_i = 0;
for i = 1:length(GLs)
    for j = 1:n_simu
        plot_md_i = plot_md_i + 1;
        if plot_md_i > 9
            nexttile(T2)
        else
            nexttile(T1)
        end
        imagesc(x,y,log10(vels{i,j})); hold on
        clim([cb_btm, cb_top])
        set(gca, 'xtick',[])
        set(gca, 'ytick',[])
        % add grounding line position
        scatter(gls{i,j}.x, gls{i,j}.y,8,'filled','r'); hold off
    end
end
colorbar_ticks = cb_btm:cb_top;
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = colorbar_ticks;
cb.TickLabels = 10.^colorbar_ticks;
cb.Label.String = 'Meter per year';
%exportgraphics(f, ['plots/steady_state/',fig_name],'ContentType','vector')
save_dir = ['plots/steady_state/',fig_name];
export_fig(save_dir)
%% Thickness
fig_name = 'h.svg';
f = figure('Position',[100,100,1000,500]);
p1=uipanel('parent',f);
p2=uipanel('parent',f);
p1.Title=("Shallow"); p1.BackgroundColor = 'w';
p2.Title=("Deep");    p2.BackgroundColor = 'w';
p1.Position=[0,0.5,1,0.5]; 
p2.Position=[0,0,1,0.5];
T1=tiledlayout(p1,3,3,'TileSpacing','none');
T2=tiledlayout(p2,3,3,'TileSpacing','none');

plot_md_i = 0;
for i = 1:length(GLs)
    for j = 1:n_simu
        plot_md_i = plot_md_i + 1;
        if plot_md_i > 9
            nexttile(T2)
        else
            nexttile(T1)
        end
        imagesc(x,y,hs{i,j}); hold on
        clim([0,800])
        set(gca, 'xtick',[])
        set(gca, 'ytick',[])
        % add grounding line position
        scatter(gls{i,j}.x, gls{i,j}.y,8,'filled','r'); hold off
    end
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Meter';
%exportgraphics(gcf, ['plots/steady_state/',fig_name],'ContentType','vector')
save_dir = ['plots/steady_state/',fig_name];
export_fig(save_dir)
%% basal drag
fig_name = 'tau_b.svg';
f = figure('Position',[100,100,1000,500]);
p1=uipanel('parent',f);
p2=uipanel('parent',f);
p1.Title=("Shallow"); p1.BackgroundColor = 'w';
p2.Title=("Deep");    p2.BackgroundColor = 'w';
p1.Position=[0,0.5,1,0.5]; 
p2.Position=[0,0,1,0.5];
T1=tiledlayout(p1,3,3,'TileSpacing','none');
T2=tiledlayout(p2,3,3,'TileSpacing','none');
% colorbar range
cb_up = 5.5; cb_btm = 3;

plot_md_i = 0;
for i = 1:length(GLs)
    for j = 1:n_simu
        plot_md_i = plot_md_i + 1;
        if plot_md_i > 9
            nexttile(T2)
        else
            nexttile(T1)
        end
        imagesc(x,y,log10(tau_bs{i,j})); hold on
        clim([cb_btm, cb_up])
        set(gca, 'xtick',[])
        set(gca, 'ytick',[])
        % add grounding line position
        scatter(gls{i,j}.x, gls{i,j}.y,8,'filled','r'); hold off
    end
end
colorbar_ticks = cb_btm:cb_up;
cb = colorbar;
cb.Ticks = colorbar_ticks;
cb.Layout.Tile = 'east';
cb.Label.String = 'Pascal';
cb.TickLabels = 10.^colorbar_ticks;
%exportgraphics(gcf, ['plots/steady_state/',fig_name],'ContentType','vector')
save_dir = ['plots/steady_state/',fig_name];
export_fig(save_dir)
%% Appendix: function
function [min_data, max_data, mean_data] = crop_domain_stats(data, W, L, front_x, ds, x, y)
%CROP_DOMAIN_STATS This function takes the given data (defined on a regular
%mesh), crops out the central portion specified by the width and length, and calculate
%the statistics of the data in the region (min, max, and mean values). x is
%the along flow direction (left to right) and y is the across-flow
%direction (bottom to top)

    % W supplied is the fjord width; the ice stream width is roughly 0.7 W
    % see .par file (parameterizing model geometry domain) for reference
    W = 0.5*W;
    % locate flow domain center line
    if rem(length(y), 2) == 0
        mid_i = length(y)/2;
    else
        mid_i = (length(y)+1)/2;
    end

    % locate the gridpoint at the front
    front_xi = find(x < front_x, 1, 'last');
    n_L = floor(L/ds);
    n_half_W = floor(W/2/ds);
    % crop
    data_c = data(mid_i - n_half_W:mid_i + n_half_W, front_xi - n_L:front_xi);

    % get stats
    mean_data = mean(data_c, 'all');
    min_data  = min(data_c, [], 'all');
    max_data  = max(data_c, [], 'all');
end