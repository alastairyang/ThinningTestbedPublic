%% Visualize the geometry and velocity field for selected glaciers at steady state
% The script recreates plots in Fig. 1
% Author: Donglai Yang
% Date: July 3, 2023

%% Parameter
ds = 50;
vel_up = 6e3; % upper limit for velocity color bar
%% Load data
load('plots/colormap/lajolla.mat');
exptname = 'MISMIP_yangTransient_Steadystate_Extended.mat';
modelnames = ["model_W5000_GL400_FC30000","model_W5000_GL400_FC120000",...
             "model_W11000_GL400_FC30000","model_W5000_GL0_FC30000"];

figure('Position',[100,100,1600,1600])
tiledlayout(2,2,"TileSpacing","none")
for md_i = 1:length(modelnames)
    mdname_char = convertStringsToChars(modelnames(md_i));
    md_dir = ['long_models_yang/' mdname_char '/' exptname];
    load(md_dir)
    
    % interp geometry to grid
    s = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, ...
                     md.results.TransientSolution(end).Surface, ds);
    base = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, ...
                     md.results.TransientSolution(end).Base, ds);
    bed = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, ...
                     md.geometry.bed, ds);

    % interp velocity to grid
    [v,x,y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, ...
                     md.results.TransientSolution(end).Vel, ds);
    [X,Y] = meshgrid(x,y);
    
    % make 3D plot
    nexttile
    bedf  = surf(X, Y, bed, 'EdgeColor','none','FaceColor','k','FaceAlpha',0.2,'EdgeAlpha',0.2); hold on
    basef = surf(X, Y, base,'EdgeColor','none','FaceColor','k','FaceAlpha',0.5); hold on
    sf = surf(X, Y, s, v,'EdgeColor','none'); hold on;
    view([130,10])
    axis off; grid off;
    clim([0, vel_up]); colormap(lajolla);
end
exportgraphics(gcf,'plots/steadystate_3Dgeom.png','Resolution',600)

%% Create a colorbar
figure('Position',[100,100,800,800])
col_data = 1e4*rand(10);
clim([0,vel_up]);
colormap(lajolla); cb = colorbar; 
cb.Label.String = "Velocity (m/a)";
cb.FontSize = 16; cb.Label.FontSize = 18;
exportgraphics(gcf,'plots/fig1_colorbar.png','Resolution',600);