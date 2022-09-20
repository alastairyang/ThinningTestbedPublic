%% Model name and parameter
model_name = 'temp_calved_md.mat';
md = load(model_name).md;

nt = size(md.results.TransientSolution,2);

% geometry parameters
Lx = max(md.mesh.x);
Ly = max(md.mesh.y);
ds = 250; % spacing, 250 meter
x = 0:ds:Lx;
y = 0:ds:Ly;
[X,~] = meshgrid(x, y);


%% diff
% Levelset and thickness
plotmodel(md, 'data',md.results.TransientSolution(end).MaskOceanLevelset - ...
                     md.results.TransientSolution(1).MaskOceanLevelset,...
              'data',md.results.TransientSolution(end).Thickness - ...
                     md.results.TransientSolution(1).Thickness);
%% just the Velocity
plotmodel(md, 'data', md.results.TransientSolution(1).Vel,...
              'data', md.results.TransientSolution(end).Vel)
%% Visualize mean velocity over time
vel_mean = zeros(nt, 1);
for i = 1:nt
    vel_mean(i) = mean(md.results.TransientSolution(i).Vel,'all');
end
figure; plot(1:nt, vel_mean)

%% Visualize ice volume and above floatation volume
vol = zeros(nt, 1);
above_vol = zeros(nt, 1);
for i = 1:nt
    vol(i) = md.results.TransientSolution(i).IceVolume;
    above_vol(i) = md.results.TransientSolution(i).IceVolumeAboveFloatation;
end
figure; 
subplot(1,2,1); plot(1:nt, vol);
subplot(1,2,2); plot(1:nt, above_vol);


%% Visualize thickness along the thalweg

if rem(size(X,1), 2) == 0
    mid_i = size(X,1)/2;
else
    mid_i = (size(X,1)+1)/2;
end
thalweg_x = X(mid_i,:);

figure;
thalweg_h_mean_line = [];
for i = 1:nt
    subplot(1,2,1)
    thalweg_h_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                          md.results.TransientSolution(i).Surface,...
                                          x, y, NaN);
    thalweg_h_line = thalweg_h_grid(mid_i,:);
    plot(thalweg_x, thalweg_h_line)
    thalweg_h_mean_line = [thalweg_h_mean_line; mean(thalweg_h_line,'all','omitnan')];
    
    pause(0.2)
end
subplot(1,2,2)
plot(1:nt, thalweg_h_mean_line)

%% Visualize thickness at sampled points along the thalweg

sample_i = 1:5:100;
figure;
thalweg_sample_ht = [];
for i = 1:nt
    thalweg_h_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                        md.results.TransientSolution(i).Surface,...
                                        x, y, NaN);
    thalweg_sample_h  = thalweg_h_grid(mid_i,sample_i);
    thalweg_sample_ht = [thalweg_sample_ht; thalweg_sample_h];
end
% Normalize the time series
thalweg_sample_ht = thalweg_sample_ht./thalweg_sample_ht(1,:);
plot(1:nt, thalweg_sample_ht)

%% Animate change in lateral profile
[geometry, ~] = query_data(index, type);
syn = testbed_data(geometry{1});
if rem(size(syn.X,1), 2) == 0
    mid_i = size(syn.X,1)/2;
else
    mid_i = (size(syn.X,1)+1)/2;
end
thalweg_x = syn.X(mid_i,:);

figure;
selected = 1:floor(nt*0.05):nt;
for i = selected
    thalweg_h_grid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x, md.mesh.y,...
                                        md.results.TransientSolution(i).Surface,...
                                        syn.x, syn.y, NaN);
    thalweg_base_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                        md.results.TransientSolution(i).Surface - md.results.TransientSolution(i).Thickness,...
                                        syn.x, syn.y, NaN);
    thalweg_h_line = thalweg_h_grid(mid_i,:);
    thalweg_base_line = thalweg_base_grid(mid_i,:);
    thalweg_both = [thalweg_h_line; thalweg_base_line];
    plot(thalweg_x, thalweg_both); hold on    
    pause(0.2)
end
%% Find grounding line
plotmodel(md, 'data', md.results.TransientSolution(1).MaskOceanLevelset <0,...
              'data', md.results.TransientSolution(end).MaskOceanLevelset <0)
          
%% Find ice presence
plotmodel(md, 'data', md.results.TransientSolution(1).MaskIceLevelset <0,...
              'data', md.results.TransientSolution(end).MaskIceLevelset <0)
%% animate groudning line
for i = 1:floor(nt*0.05):nt
    plotmodel(md, 'data', md.results.TransientSolution(i).MaskOceanLevelset <0)
    pause(0.1)
end

%% Animate the Surface Elevation
for i = 1:floor(nt*0.05):nt
    plotmodel(md, 'data', md.results.TransientSolution(i).Surface, 'caxis', [0,800])
    pause(0.1)
end
%% animate velocity
for i = 1:floor(0.05*nt):nt
    plotmodel(md, 'data',md.results.TransientSolution(i).Vel);
    pause(0.1)
end
%% Animate change in surface elevation
for i = 1:nt-1
    plotmodel(md, 'data', md.results.TransientSolution(i+1).Surface - md.results.TransientSolution(i).Surface, 'caixs', [-30, 5])
    pause(0.1)
end

%% Apply EOF to the surface signals
% We collect all data into a 3-D matrix
% Before that, we need to re-grid into meshgrid
x_max = max(md.mesh.x);
x_min = min(md.mesh.x);
y_max = max(md.mesh.y);
y_min = min(md.mesh.y);
ns = 150; % same resolution as BedMachine; not necessary tho
x = x_min:ns:x_max;
y = y_min:ns:y_max;
[X, Y] = meshgrid(x,y);

Ht = zeros(numel(y), numel(x), nt);
for i = 1:nt
    this_Ht = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x, md.mesh.y,...
                       md.results.TransientSolution(i).Surface,...
                       x, y);
    Ht(:,:,i) = this_Ht;
end

% EOF
[eofmaps, pc, expvar] = eof(Ht,5);

% Plot
figure; subplot(2,2,1); imagesc(x,y,eofmaps(:,:,1)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 1')
subplot(2,2,2); imagesc(x, y,eofmaps(:,:,2)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 2')
% add the grounding line and bed plot
bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                       md.geometry.bed,...
                       x, y);
levelset = InterpFromMeshToGrid(md.mesh.elements. md.mesh.x, md.mesh.y,...
                       md.results.TransientSolution(end).MaskOceanLevelset,...
                       x, y);
subplot(2,2,3); imagesc(x, y, levelset); axis xy image off; title('Ocean Levelset')
set(gca, 'CLim', [-5,5])
subplot(2,2,4); imagesc(x, y, bed);  axis xy image off; title('Bed')

%% Plot reconstructed time series from principal components
t = 0:0.1:(nt-1)*0.1;
figure;
subsubplot(5,1,1)
plot(t, pc(1,:)/pc(1,1))
box off
axis tight
ylabel(['pc1: ',num2str(expvar(1)),'%'])
title 'The first two PCs'

subsubplot(5,1,2)
plot(t, pc(2,:)/pc(2,1))
box off
axis tight
set(gca,'yaxislocation','right')
ylabel(['pc2: ',num2str(expvar(2)),'%'])

subsubplot(5,1,3)
plot(t, pc(3,:)/pc(3,1))
box off
axis tight
ylabel(['pc3: ',num2str(expvar(3)),'%'])

subsubplot(5,1,4)
plot(t, pc(4,:)/pc(4,1))
box off
axis tight
set(gca,'yaxislocation','right')
ylabel(['pc4: ',num2str(expvar(4)),'%'])

subsubplot(5,1,5)
plot(t, pc(5,:)/pc(5,1))
box off
axis tight
ylabel(['pc5: ',num2str(expvar(5)),'%'])

figure;
% reconstruct one time series from the first 5 PCs
sum_h = zeros(size(t));
for i = 1:5
    sum_h = sum_h + pc(i,:)*0.01*expvar(i);
end
plot(t, sum_h)
%% plot time series at all points
coor_N = numel(x)*numel(y);
t = 0:0.1:(nt-1)*0.1;
all_timeseries = zeros(coor_N, nt);
figure
for i = 1:coor_N
    % get the corresponding subscript index
    [row, col] = ind2sub([numel(y), numel(x)], i);
    this_timeserie = Ht(row, col,:);
    this_timeserie = this_timeserie(:);
    this_timeserie = this_timeserie/this_timeserie(1);
    all_timeseries(i,:) = this_timeserie;
    plot(t, this_timeserie)
    hold on;
end

    