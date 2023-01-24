%% Map view: analyze localized basal perturbation
% here we use mode decomposition to decompose the H(x,y,t) from the
% localized basal perturbation experiment. The goal is to arrive at a
% simple clean data reduction that shows the impact of the localized basal
% perturbation: the magnitude of cyclic and trend change, the propagation
% across the glacier, and etc

%% Experiment with point-wise timeseries analysis
md_expt = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat').md;
md_ctrl = load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
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
md_grid_v = reshape(md_grid, [xl*yl,size(md_grid,3)]);
STs = detrend(md_grid_v, 6);
LTs = md_grid_v - STs;
LTs = reshape(LTs, [xl, yl, nt]);
STs = reshape(STs, [xl, yl, nt]);

% visualize
figure;for i = 1:nt; imagesc(x,y,squeeze(LTs(:,:,i)));clim([-10,0]);pause(0.05);end
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