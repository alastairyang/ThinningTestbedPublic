%% Plot the map view 2D velocity, thickness, basal drag at the steady state
geom_type = 'deep'; % options: "deep" or "shallow"
ds = 50; % grid size for the regular grid

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
% iterate over Deep GL models
%tiledlayout(3,3,'TileSpacing','tight')
switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

n_simu = size(folder_dir_groups{geom_i}, 1);
vels = cell(1, n_simu);
hs = cell(1, n_simu);
tau_bs = cell(1, n_simu);
gls = cell(1, n_simu);

for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md = load([group.folder{j},'/', group.name{j}, '/', md_name]).md;
    % query data from the last time step
    vel = md.results.TransientSolution(end).Vel;
    h   = md.results.TransientSolution(end).Thickness;
    % get the sliding law coefficient and calculate the basal drag
    tau_b = md.friction.C.^2.*(vel/md.constants.yts);
    % interp to regular grid
    [vels{j}, ~, ~]  = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, vel, ds);
    [hs{j}, ~, ~]    = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, h, ds);
    [tau_bs{j}, ~, ~]= mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, tau_b, ds);
    % get the grounding line positions
    gls{j} = isoline(md, md.results.TransientSolution(end).MaskOceanLevelset,'value',0);
    % if multiple enclosed areas were found, we only save the one that has
    % most nodes
    if size(gls{j},2) ~= 1
        gl_tbl = struct2table(gls{j});
        gl_tbl = sortrows(gl_tbl, 'nods','descend');
        gls{j} = table2struct(gl_tbl(1,:));
    end
        
end
% get the x,y axis
[~, x, y]= mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, tau_b, ds);

% make a tiled plot
%% velocity
fig_name = ['vel_',geom_type,'.pdf'];
figure('Position',[100,100,1000,500])
tiledlayout(3,3, "TileSpacing","none")
for j = 1:n_simu
    nexttile
    imagesc(x,y,log10(vels{j})); hold on
    clim([0,4])
    set(gca, 'xtick',[])
    set(gca, 'ytick',[])
    % add grounding line position
    scatter(gls{j}.x, gls{j}.y,8,'filled','r'); hold off
end
colorbar_ticks = 0:4;
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = colorbar_ticks;
cb.TickLabels = 10.^colorbar_ticks;
exportgraphics(gcf, ['plots/steady_state/',fig_name],'ContentType','vector')

%% Thickness
fig_name = ['h_',geom_type,'.pdf'];
figure('Position',[100,100,1000,500])
tiledlayout(3,3, "TileSpacing","none")
for j = 1:n_simu
    nexttile
    imagesc(x,y,hs{j}); hold on
    clim([0,800])
    set(gca, 'xtick',[])
    set(gca, 'ytick',[])
    % add grounding line position
    scatter(gls{j}.x, gls{j}.y,8,'filled','r'); hold off
end
cb = colorbar;
cb.Layout.Tile = 'east';
exportgraphics(gcf, ['plots/steady_state/',fig_name],'ContentType','vector')

%% basal drag
fig_name = ['taub_',geom_type,'.pdf'];
figure('Position',[100,100,1000,500])
tiledlayout(3,3, "TileSpacing","none")
for j = 1:n_simu
    nexttile
    imagesc(x,y,tau_bs{j}); hold on
    clim([5e4,3e5])
    set(gca, 'xtick',[])
    set(gca, 'ytick',[])
    % add grounding line position
    scatter(gls{j}.x, gls{j}.y,8,'filled','r'); hold off
end
cb = colorbar;
cb.Layout.Tile = 'east';
exportgraphics(gcf, ['plots/steady_state/',fig_name],'ContentType','vector')
