%% Map view: analyze localized basal perturbation
% here we use mode decomposition to decompose the H(x,y,t) from the
% localized basal perturbation experiment. The goal is to arrive at a
% simple clean data reduction that shows the impact of the localized basal
% perturbation: the magnitude of cyclic and trend change, the propagation
% across the glacier, and etc

%% Read in one model
load('long_models_yang/model_W11000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat');

% we need to transform to regular grid
[md_grid, x, y] = mesh_to_grid_overtime(md.mesh.elements, md.mesh.x, md.mesh.y, results_tbl.Thickness, 50);

% DMD
[evals, evecs, mode_amps, mode_freq, growth_rates, mode_E] = dmd_rom(md_grid, 5, 0.1);