%% Analysis of effective pressure experiment at a wavy/rough bed
% This script creates figure 5 in the main text
% Author: Donglai Yang
% Date: June 27, 2023
clear; clc;
% Load model data, "md"
load('wavybed_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat')

%% Extract data
sample_interval = 400; % meter; distance between two control points, if there are multiple sampled points
sample_number = 1; % make # samples along the center line
ds = 50; % grid spacing, meter
ds_i = sample_interval/ds; % spacing in grid index

% # of timesteps
nt = size(md.results.TransientSolution,2);

% save grounding line and front locations
gl_locs = zeros(1,nt);
front_locs = zeros(1,nt);
for k = 1:nt
    front_mask = md.results.TransientSolution(k).MaskIceLevelset;
    gl_mask = md.results.TransientSolution(k).MaskOceanLevelset;
    front_locs(k) = locate_calvingfront(md, front_mask);
    gl_locs(k) = locate_groundingline(md, gl_mask);
end

% Get the thickness, control point location, and time axis
% Create regular meshgrid. We will interp from unstructured mesh to regular grid
% later
Lx = max(md.mesh.x);
Ly = max(md.mesh.y);
x = 0:ds:Lx-ds;
y = 0:ds:Ly-ds;
[X,Y] = meshgrid(x,y);
if rem(size(X,1), 2) == 0
    mid_i = size(X,1)/2-10;
else
    mid_i = (size(X,1)+1)/2-10;
end
thalweg_x = X(mid_i,:);
mesh_elements = md.mesh.elements;
mesh_x = md.mesh.x;
mesh_y = md.mesh.y;
% retrieve the last position of the ice front; spclevelset ~= 0
pos = find(md.levelset.spclevelset(1:end-1,end) < 1 & md.levelset.spclevelset(1:end-1,end) > -1);
x_front = min(mesh_x(pos));

% find grid index where we want to sample h(t)
[~, x_i_nearest] = min(abs(x - x_front));
% sampled points: start at 1 spacing distance behind the last ice front
front_i = x_i_nearest-ds_i;
end_i   = front_i - sample_number*ds_i;
sample_i = front_i:-ds_i:(end_i+ds_i);
% get the absolute distance
sample_x = x(sample_i);
sample_y = max(y) - ds*mid_i;

% Convert structure table to table, for easier data wrangling
calve_results = struct2table(md.results.TransientSolution);
empty_md = md;
empty_md.results.TransientSolution = [];

% initialize space to store h(t)
thalweg_sample_ht = [];
thalweg_sample_st = [];

% get all h(t) and s(t) data
for j = 1:nt
    % from mesh to grid
    calve_surface = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
        calve_results.Surface{j},x, y, NaN);
    calve_thickness = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
        calve_results.Thickness{j},x, y, NaN);
    surface_profile  = calve_surface(mid_i,sample_i);
    thickness_profile = calve_thickness(mid_i,sample_i);
    thalweg_sample_ht = [thalweg_sample_ht; thickness_profile];
    thalweg_sample_st = [thalweg_sample_st; surface_profile];
end
thalweg_sample_ht = thalweg_sample_ht - thalweg_sample_ht(1,:);
thalweg_sample_st = thalweg_sample_st - thalweg_sample_st(1,:);
% crop the one extra time step in calving simulation
thalweg_sample_ht = thalweg_sample_ht(1:end-1,:);
thalweg_sample_st = thalweg_sample_st(1:end-1,:);

% time vector
time = calve_results.time;
% shift time to start at 0
time = time(1:end-1);
time = time - time(1);

% save data
ht_data.h = thalweg_sample_ht;
ht_data.s = thalweg_sample_st;
ht_data.t = time;
ht_data.x = sample_x;

htt = diff(thalweg_sample_ht,1,1);
stt = diff(thalweg_sample_st,1,1);

% use the following command to plot all h(t) and dh/dt(t) timeseries
% figure; plot(htt) 
% figure; plot(thalweg_sample_ht)

% Load wavy bed data
wavybed = load('random_beds/dD80_H07_RO4000_3_grid.mat').wavybed;

% plot the glacier lateral profile evolution
% shift data to the start of perturbation
idx = time>5;
plot_t = time(idx)-5;

thalweg_sample_ht_c = thalweg_sample_ht(idx);
front_locs = front_locs(1:end-1); front_locs = front_locs(idx);
gl_locs = gl_locs(1:end-1); gl_locs = gl_locs(idx);

% load colormap
cp = load('plots/colormap/lajolla.mat').lajolla;
cp_modi = [cp'*255; ones(1,size(cp,1))];
cp_modi_interp = transpose(uint8(interp1(1:256, cp_modi', (size(cp,1)/length(plot_t))*(1:length(plot_t)))));
brown = [228, 119, 41]/255;

% highlight times of stepwise retreat with circular markers
t_start = 6.5;
t_end = 12.5;
t_mark = t_start:0.05:t_end;
gl_mark = interp1(plot_t, gl_locs/1000, t_mark);

%% make plots
figure('Position',[100,100,800,500])
plot_t_max = 16; % only plot the first 16 years (full perturbation; no relaxation plotted)
plot_idx = plot_t<plot_t_max;
% thickness change timeseries
tiledlayout(2,1,'TileSpacing','tight')
nexttile
yyaxis left;
p1 = plot(plot_t(plot_idx), thalweg_sample_ht_c(plot_idx),'-k','LineWidth',1.5);hold on;
%colororder(cool(sample_number)); 
ylabel('Meter','FontSize',13)
set(gca,'ycolor','k')
hold on;

yyaxis right;
% plot front and grounding line location
p2 = plot(plot_t(plot_idx), front_locs(plot_idx)/1000,'-.','LineWidth',1.5,'Color',brown);hold on;
p3 = plot(plot_t(plot_idx), gl_locs(plot_idx)/1000,'-.k','LineWidth',1.5);
hold on; 
%legend(["Thickness","Front","Grounding Line"],'FontSize',13,'Location','southwest')
%legend([p1,p2,p3],'Thickness','Front','Grounding Line','FontSize',13,'Location','southwest')
set(gca,'ycolor',brown)
ylabel('Kilometer','FontSize',13)
xlabel('Year since perturbation starts','FontSize',13)
xlim([0,plot_t_max])

% plot profile evolution
nexttile
[x_plot, bed_plot,~,~,colors_plot,t] = plot_all_profiles(md,10,cp); hold on;
% add a dashed line at where the timeseries is taken
xline(sample_x/1000,':k','LineWidth',2)
xlabel('Along flow direction (km)','FontSize',13);
ylabel('Elevation (m)','FontSize',13);
% add markers
c_markers = interp1(t, colors_plot, t_mark);
gl_at_base = interp1(x_plot/1e3, bed_plot, gl_mark);
scatter(gl_mark, gl_at_base, 70, c_markers,"filled",'MarkerEdgeColor','k','MarkerEdgeAlpha',0.02); 

% back to the first subplot and add markers
nexttile(1)
yyaxis right
scatter(t_mark, gl_mark, 70, c_markers,'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.02)
legend([p1,p2,p3],'Thickness','Front','Grounding Line','FontSize',13,'Location','southwest')

exportgraphics(gcf, 'plots/wavybed_gl_retreat_line.png','Resolution',600)

% Create a colorbar for the profiles
val_max = 10;
fake_data = 26*rand(10,10);
figure; imagesc(fake_data); 
clim([0,val_max]); 
colormap(colors_plot);
cb_p = colorbar;
cb_p.Ticks = [0,floor(val_max/2),val_max];
cb_p.TickLabels = string(cb_p.Ticks);
cb_p.FontSize = 14;
cb_p.Location = 'northoutside';
%exportgraphics(gcf,'plots/wavy_colorbar.pdf','ContentType','vector')
%% Map view of bed topography and thinning rate
% crop the plan view extent (zoom in)
% start by specifying x and y limits of the cropped data
x_lim = [3e4, 5e4];
y_lim = [3e3, 9e3];
xi_keep = x > x_lim(1) & x < x_lim(2);
yi_keep = y > y_lim(1) & y < y_lim(2);
x_c = x(xi_keep);
y_c = y(yi_keep);
[X_c, Y_c] = meshgrid(x_c, y_c);

% Plan view
i = 260; % the index of the time slice (last time slice)
dh = md.results.TransientSolution(i).Surface - md.results.TransientSolution(i-1).Surface;
dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
dhdt = dh/dt;
time = md.results.TransientSolution(i-1).time;
% load colormap
load('plots/colormap/nuuk_polar.mat')
[dhdt,~,~] = mesh_to_grid(md.mesh.elements,md.mesh.x,md.mesh.y,dhdt,ds);

fig = figure('Position',[100,100,800,250]);
% bed topography as image ('imagesc'); thinning rates as contour
% ('contour')
ax1 = axes(fig); 
ax2 = copyobj(ax1,fig);
contour(ax1, X_c/1000, Y_c/1000, dhdt(yi_keep, xi_keep),"ShowText",true,"LabelFormat","%0.0f m/a",'LineWidth',1.3,'LineStyle','-');
clim([-12,12]);colormap(ax1, nuuk);
xlabel('Along flow direction (km)','FontSize',13);hold on;
imagesc(ax2, x_c/1000, y_c/1000, wavybed(yi_keep, xi_keep),'AlphaData',0.6)
colormap(ax2, nuuk);clim(ax2,[-450,-200]);cb = colorbar(ax2);
cb.Label.String = 'Depth (m)';
cb.Label.FontSize = 13;
hold on

% add red point to where the timeseries is taken
scatter(ax1,sample_x/1000, sample_y/1000,400,'r','filled');

ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';

yticks([2,4,6,8])
xticks([35,40,45])
exportgraphics(gcf, 'plots/wavybed_gl_retreat_map.png','Resolution',600)

%% functions
function [x,bed_profile, surface_profiles, base_profiles, colors_p, t_sample] = plot_all_profiles(md, skip, colormap_p)
%PLOT_ALL_PROFILES
%   Plot the lateral profiles of the glaciers at specified time steps
%
%   Input:
%       md        [ISSM model]: ISSM model class
%       skip      [int]: the number of timesteps we skip in processing
%       colormap_p[array]: Colormap array, usually 256x3
%       
%   Output:
%       x [array]: 1d vector of x axis
%       bed_profiles     [array]: bed profile
%       surface_profiles [array]: surface elevation at sampled timesteps
%       base_profiles    [array]: base elevation at sampled timesteps

    if nargin == 1; skip = 1; end
    nt = size(md.results.TransientSolution,2);
    t = [md.results.TransientSolution(:).time];
    t = t-t(1);
    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)

    % color gradient
    if nargin < 3 % if colormap not provided
        color_length = nt;
        red = [255, 51, 153]/255;
        sth = [153, 153, 255]/255;
        colors_p = [linspace(red(1),sth(1),color_length)',...
                    linspace(red(2),sth(2),color_length)',...
                    linspace(red(3),sth(3),color_length)'];
    else
        % create color order from colormap
        n_c = length(1:skip:nt);
        colors_p = colormap_to_colororder(colormap_p, n_c,1,100);
        colors_p2 = repelem(colors_p, 2*ones(n_c,1), 1);
    end
    % bedrock color
    rock_rgb = [204, 204, 204]/255;
    % y base
    basevalue = -600;
    
    % geometry parameters
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 250; % spacing, 250 meter
    x = 0:ds:Lx-ds;
    y = 0:ds:Ly-ds;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);

    % preallocate
    surface_profiles = [];
    base_profiles = [];

    % sample some time slices
    t_sample_idx = 1:skip:nt;
    t_sample = t(t_sample_idx);

    % plot bed profile first
    for i = t_sample_idx
        % mesh to grid
        surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Surface,...
            x, y, NaN);
        base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Base,...
            x, y, NaN);
        surface_profile  = surface(mid_i,:);
        base_profile = base(mid_i,:);
        surface_profiles = [surface_profiles; surface_profile];
        base_profiles = [base_profiles; base_profile];
        % plot
        plot(x/1000, surface_profile, 'LineWidth',1.2); hold on
        plot(x/1000, base_profile, 'LineWidth',1.2);hold on
        
        ylim([basevalue,1500])
    end
    plot(x/1000, bed_profile, 'black');hold on;
    colororder(colors_p2)
    % color the bedrock
    a = area(x/1000, bed_profile, basevalue); hold off
    a.FaceColor = rock_rgb; a.EdgeColor = rock_rgb;

end