sample_interval = 400; % meter; distance between two control points
sample_number = 1; % make # samples along the center line
ds = 50; % grid spacing, meter
ds_i = sample_interval/ds;

%md = load([folder_dir(i).folder '/' fullname]).md;
nt = size(md.results.TransientSolution,2);

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
x = 0:ds:Lx;
y = 0:ds:Ly;
[X,~] = meshgrid(x, y);
if rem(size(X,1), 2) == 0
    mid_i = size(X,1)/2-10;
else
    mid_i = (size(X,1)+1)/2-10;
end
thalweg_x = X(mid_i,:);
% save the mesh elements and (x,y)
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

%% plot
figure('Position',[100,100,600,600])
yyaxis left; 
plot(ht_data.t(1:end-1), htt,':k','LineWidth',1.5);colororder(cool(sample_number)); hold on;
plot(ht_data.t(1:end-1), stt,'-k','LineWidth',1.5);colororder(cool(sample_number)); hold on;
hold on; 
yyaxis right;plot(ht_data.t(1:end-1), front_locs(1:258),'-.r','LineWidth',1.5);hold on;plot(ht_data.t(1:end-1),gl_locs(1:258),'-.b','LineWidth',1.5)
legend(["dh/dt","ds/dt","Front","GL"],'FontSize',13,'Location','southwest')
exportgraphics(gcf,'plots/wavybed_gl_front_dhdt.png','Resolution',300)
%%
i = 200;
dh = md.results.TransientSolution(i).Surface - md.results.TransientSolution(i-1).Surface;
dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
dhdt = dh/dt;
time = md.results.TransientSolution(i-1).time;
load('plots/colormap/lajolla.mat')
[dhdt,~,~] = mesh_to_grid(md.mesh.elements,md.mesh.x,md.mesh.y,dhdt,ds);

figure('Position',[100,100,1000,150]);imagesc(x,y,dhdt);clim([-20,0]);colormap(lajolla);colorbar

