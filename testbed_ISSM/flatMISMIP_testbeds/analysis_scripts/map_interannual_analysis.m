%% Map view: analyze localized basal perturbation
% here we use mode decomposition to decompose the H(x,y,t) from the
% localized basal perturbation experiment. The goal is to arrive at a
% simple clean data reduction that shows the impact of the localized basal
% perturbation: the magnitude of cyclic and trend change, the propagation
% across the glacier, and etc

%% Experiment with point-wise timeseries analysis
% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat';
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);
% plot parameter 
ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

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

% shallower grounding line
n_simu = size(folder_dir_groups{shallowGL_i}, 1);
% pre-allocate
deltaH_ctrl = cell(n_simu, 2);
deltaH_expt = cell(n_simu, 2);
gl_cells_ctrl = cell(n_simu,2);
gl_cells_expt = cell(n_simu,2);

% iterate over Deep GL models
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{deeperGL_i};
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % isolate the effect from localized basal perturbation
    deltaH = [results_tbl_expt.Thickness{:}] - [results_tbl_ctrl.Thickness{:}];
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
    % timeseries decomposition
    % we reshape x,y into a long vector and after decomposed return to a map
    xl = size(md_grid,1); yl = size(md_grid,2); nt = size(md_grid,3);
    md_grid_v = transpose(reshape(md_grid, [xl*yl, size(md_grid,3)]));
    STs = detrend(md_grid_v, 5);
    LTs = md_grid_v - STs;
    LTs = reshape(LTs, [nt, xl, yl]);
    STs = reshape(STs, [nt, xl, yl]);
    % permute back (time axis at the 3rd dimension)
    LTs = permute(LTs, [2,3,1]);
    STs = permute(STs, [2,3,1]);
    
    % animate
%     figure;for i = 1:nt; imagesc(x,y,squeeze(LTs(:,:,i)));clim([-10,0]);pause(0.05);end
%     figure;for i = 1:nt; imagesc(x,y,squeeze(STs(:,:,i)));clim([-10,0]);pause(0.05);end
    
    % find the spatial mean and spatial std
    LTs_crop = LTs(:,:,51:end); STs_crop = STs(:,:,51:end); % skip the first 5 years
    LTs_mean = mean(LTs_crop, 3); STs_max = max(STs_crop,[], 3);STs_min = min(STs_crop,[], 3);
    LTs_std  = std(LTs_crop, 0, 3);  STs_std  = std(STs_crop, 0, 3);
    % long term mean timeseries
    LTs_time_mean = squeeze(nanmean(LTs, [1,2]));
    
    % Visualization
    % long term: Total thinning extent and mean timseries 
    figure('Position',[100,100,1100,400]);
    t = tiledlayout(2,2);
    md_name = md_ctrl.miscellaneous.name(9:end);
    title(t,md_name,'interpreter','none')
    nexttile;
    imagesc(x,y,squeeze(LTs(:,:,end))); colorbar; clim([-15,0]);
    title('Total thinning from trend component'); ylabel('Meter')
    nexttile;
    plot(0:0.1:26-0.1, squeeze(LTs_time_mean)); xlim([0,26])
    title('Thinning trend'); xlabel('Year'); 
    nexttile;
    imagesc(x,y,STs_max - STs_min);colorbar;clim([-5,5]);
    title('Mean cyclic magnitude'); 
    nexttile;
    imagesc(x,y,STs_std); colorbar;clim([0,5]);
    title('1 std cyclic mangitude')

    % save the plot
    save_dir = ['plots/localized_basal_perturb/',md_name,'.png'];
    exportgraphics(gcf, save_dir, 'Resolution',300)
    disp(['model ',md_name,' is completed!'])
end

%% Experiment with EOF
% parameter
n_m = 5;

md_expt = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat').md;
md_ctrl = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
deltaH = [results_tbl_expt.Thickness{:}] - [results_tbl_ctrl.Thickness{:}];
deltaH_cell = num2cell(deltaH,1);

% h_1st = results_tbl.Thickness(2:end);
% h_2nd = results_tbl.Thickness(1:end-1);
% dhdt = [h_2nd{:}] - [h_1st{:}];
% % wrap into cell array
% dhdt_cell = num2cell(dhdt,1);
% clear dhdt
% we need to transform to regular grid
[md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
% suppress large dhdt due to frontal ablation
[mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
% dh/dt
md_grid = md_grid(2:end,:,:) - md_grid(1:end-1,:,:);
mask_grid = mask_grid(2:end,:,:);
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
subplot(2,3,1); imagesc(x,y,eofmaps(:,:,1)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 1')
subplot(2,3,2); imagesc(x, y,eofmaps(:,:,2)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 2')
subplot(2,3,3); imagesc(x, y,eofmaps(:,:,3)); axis xy image off; cmocean('delta','pivot',0);
title('Mode 3')
% temporal evolution
subplot(2,3,4); anomaly(0:0.1:25.9-0.1,pc(1,:));
subplot(2,3,5); anomaly(0:0.1:25.9-0.1,pc(2,:));
subplot(2,3,6); anomaly(0:0.1:25.9-0.1,pc(3,:));

%% Experiment with DMD
% use two models as an illustration (localized basal perturbation)
% parameter
n_m = 5;

md_expt = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat').md;
md_ctrl = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
results_tbl_expt = struct2table(md_expt.results.TransientSolution);
results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
deltaH = [results_tbl_expt.Thickness{:}] - [results_tbl_ctrl.Thickness{:}];
deltaH_cell = num2cell(deltaH,1);

% h_1st = results_tbl.Thickness(2:end);
% h_2nd = results_tbl.Thickness(1:end-1);
% dhdt = [h_2nd{:}] - [h_1st{:}];
% % wrap into cell array
% dhdt_cell = num2cell(dhdt,1);
% clear dhdt
% we need to transform to regular grid
[md_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, 50);
% suppress large dhdt due to frontal ablation
[mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, 50);
md_grid(mask_grid>=0) = 0;

% DMD
[evals, evecs, mode_amps, mode_freq, growth_rates, mode_E] = dmd_rom(md_grid2, n_m, 0.1); % first 5 mode
dominant_mode_No=find(mode_amps==max(mode_amps));
%==========Mode Amplitudes===========
figure('color', 'w');
stem(abs(mode_freq),abs(mode_amps), 'k^', 'filled'); hold on
plot(abs(mode_freq(dominant_mode_No)),abs(mode_amps(dominant_mode_No)),'ro','MarkerSize',10);
set(gca, 'YScale', 'log')
title('Mode amplitude versus Frequency');
xlabel('f [Hz]');ylabel('Amplitude');
%==========Evolve one mode in time============
m=2; % plot the first mode
theta=linspace(0,6*pi,100);
figure('color', 'w','Position',[100,100,1000,300]); 
for i=1:length(theta)
    H=squeeze(real(evecs(m,:,:)*exp(1i*theta(i))));
    imagesc(x,y,H)
    clim([-0.01, 0.01])
    colorbar
    pause(0.05);
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