sample_interval = 400; % meter; distance between two control points
sample_number = 1; % make # samples along the center line
ds = 50; % grid spacing, meter
ds_i = sample_interval/ds;

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
% geometry parameters
Lx = max(md.mesh.x);
Ly = max(md.mesh.y);
x = 0:ds:Lx-ds;
y = 0:ds:Ly-ds;
[X,Y] = meshgrid(x,y);
[X,~] = meshgrid(x, y);
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
% sampled points: start at __ behind the last ice front
front_i = x_i_nearest-ds_i;
end_i   = front_i - sample_number*ds_i;
sample_i = front_i:-ds_i:(end_i+ds_i);
% get the absolute distance
sample_x = x(sample_i);
sample_y = max(y) - ds*mid_i;

% remove model class; data store in table instead to clear space
calve_results = struct2table(md.results.TransientSolution);
empty_md = md;
empty_md.results.TransientSolution = [];

% initialize space to store h(t)
thalweg_sample_ht = [];
thalweg_sample_st = [];

% get all h(t) data
for j = 1:nt
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

%% Load wavy bed data
wavybed = load('random_beds/dD80_H07_RO4000_3_grid.mat').wavybed;

%% plot
figure('Position',[100,100,800,500])
tiledlayout(2,1,'TileSpacing','tight')
%subplot(3,1,1)
nexttile
yyaxis left;
plot(ht_data.t(1:end), thalweg_sample_ht,'-b','LineWidth',1.5);colororder(cool(sample_number)); hold on;
ylabel('m','FontSize',13)
hold on;
%plot(ht_data.t(1:end-1), htt,'-k','LineWidth',1.5);colororder(cool(sample_number)); hold on;
yyaxis right;
plot(ht_data.t(1:end-1), front_locs(1:258)/1000,'-.r','LineWidth',1.5);hold on;
plot(ht_data.t(1:end-1),gl_locs(1:258)/1000,'-.b','LineWidth',1.5)
legend(["H(t)","Front","Grounding Line"],'FontSize',13,'Location','southwest')
ylabel('km','FontSize',13)
xlabel('Year','FontSize',13)
xlim([0,26])

% profile
%subplot(3,1,2)
nexttile
plot_all_profiles(md,10)
% add a dashed line at where the timeseries is taken
xline(sample_x/1000,'-.b','LineWidth',0.5)
xlabel('Along flow direction (km)','FontSize',13);
ylabel('Elevation (m)','FontSize',13);

exportgraphics(gcf, 'plots/wavybed_gl_retreat_line.png','Resolution',600)

%% Map view
% crop the plan view extent (zoom in)
x_lim = [3e4, 5e4];
y_lim = [3e3, 9e3];
xi_keep = x > x_lim(1) & x < x_lim(2);
yi_keep = y > y_lim(1) & y < y_lim(2);
x_c = x(xi_keep);
y_c = y(yi_keep);
[X_c, Y_c] = meshgrid(x_c, y_c);

% Plan view
%subplot(3,1,3)
i = 260; % the index of the time slice
dh = md.results.TransientSolution(i).Surface - md.results.TransientSolution(i-1).Surface;
dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
dhdt = dh/dt;
time = md.results.TransientSolution(i-1).time;
load('plots/colormap/nuuk_polar.mat')
[dhdt,~,~] = mesh_to_grid(md.mesh.elements,md.mesh.x,md.mesh.y,dhdt,ds);

fig = figure('Position',[100,100,800,250]);

ax1 = axes(fig); 
ax2 = copyobj(ax1,fig);
contour(ax1, X_c/1000, Y_c/1000, dhdt(yi_keep, xi_keep),"ShowText",true,"LabelFormat","%0.0f m/a",'LineWidth',1.3,'LineStyle','-');
clim([-12,12]);colormap(ax1, nuuk);
%exportgraphics(gcf,'plots/wavybed_gl_front_dhdt.png','Resolution',300)
xlabel('Along flow direction (km)','FontSize',13)
% add red point to where the timeseries is taken
hold on
imagesc(ax2, x_c/1000, y_c/1000, wavybed(yi_keep, xi_keep),'AlphaData',0.6)
colormap(ax2, nuuk);clim(ax2,[-450,-200]);cb = colorbar(ax2);
cb.Label.String = 'Depth (m)';
cb.Label.FontSize = 13;
hold on
scatter(ax1,sample_x/1000, sample_y/1000,400,'r','filled');

ax2.UserData = linkprop([ax1,ax2],...
    {'Position','InnerPosition','xtick','ytick', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
ax2.Visible = 'off';

yticks([2,4,6,8])
xticks([35,40,45])
exportgraphics(gcf, 'plots/wavybed_gl_retreat_map.png','Resolution',600)

%% functions
function plot_all_profiles(md, skip)

    if nargin == 1; skip = 1; end
    nt = size(md.results.TransientSolution,2);
    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)

    % color gradient
    color_length = nt;
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
                linspace(red(2),sth(2),color_length)',...
                linspace(red(3),sth(3),color_length)'];
    
    % geometry parameters
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 250; % spacing, 250 meter
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);

    bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);
    plot(x/1000, bed_profile, 'black');hold on;
    for i = 1:skip:nt
        
        surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Surface,...
            x, y, NaN);
        base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(i).Base,...
            x, y, NaN);
        bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
        surface_profile  = surface(mid_i,:);
        base_profile = base(mid_i,:);
        
        plot(x/1000, bed_profile, 'k');hold on;
        plot(x/1000, surface_profile, 'Color',colors_p(nt+1-i,:)); hold on
        plot(x/1000, base_profile,'Color',colors_p(nt+1-i,:));hold on
        
        %hold off
        ylim([-800,1500])
    end
    hold off
end