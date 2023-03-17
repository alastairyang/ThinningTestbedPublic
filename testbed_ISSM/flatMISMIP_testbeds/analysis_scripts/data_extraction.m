%% Analaysis script
% 1. Here we are going to export the TransientSolution into data format that
%    is easier to work with, meanwhile extracting the timeseries at selected
%    points along the centerline.

% This file should be run at "flatMISMIP_testbeds"

%% global parameters
% parameters
sample_interval = 1000; % meter; distance between two control points
sample_number = 20; % make # samples along the center line
ds = 50; % grid spacing, meter
foldername = 'long_models_yang';

% dictionary
dict.calving.modelname = 'MISMIP_yangTransient_CalvingOnly.mat';
dict.calving.save_foldername = 'analyzed_data/calve_only/';
dict.calving.save_fileprefix = 'ht_calve_';

% dict.calving_smw.modelname = 'MISMIP_yangTransient_Calving_SMweakening.mat';
% dict.calving_smw.save_foldername = 'analyzed_data/smw_calve/';
% dict.calving_smw.save_fileprefix = 'ht_smw_calve_';

dict.calving_mu.modelname = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
dict.calving_mu.save_foldername = 'analyzed_data/mu_calve/';
dict.calving_mu.save_fileprefix = 'ht_mu_calve_';

% Mass unloading + localized basal perturbation: transient pulse
dict.calving_mu_pulse8.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_PulseGaussianPerturb_8.mat';
dict.calving_mu_pulse8.save_foldername = 'analyzed_data/mu_pulse_calve/';
dict.calving_mu_pulse8.save_fileprefix = 'ht_mu_pulse_calve_';

% Mass unloading + localized basal perturbation: diffused pulse
dict.calving_mu_diffu8.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_DiffuGaussianPerturb_8.mat';
dict.calving_mu_diffu8.save_foldername = 'analyzed_data/mu_diffu_calve/';
dict.calving_mu_diffu8.save_fileprefix = 'ht_mu_diffu_calve_';

% localized basal perturbation: transient pulse
dict.calving_pulse8.modelname = 'MISMIP_yangTransient_Calving_PulseGaussianPerturb_8.mat';
dict.calving_pulse8.save_foldername = 'analyzed_data/pulse_calve/';
dict.calving_pulse8.save_fileprefix = 'ht_pulse_calve_';

% localized basal perturbation: diffused pulse
dict.calving_diffu8.modelname = 'MISMIP_yangTransient_Calving_DiffuGaussianPerturb_8.mat';
dict.calving_diffu8.save_foldername = 'analyzed_data/diffu_calve/';
dict.calving_diffu8.save_fileprefix = 'ht_diffu_calve_';

% Old experiment; localized perturbation
% dict.calving_gp1.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_1.mat';
% dict.calving_gp1.save_foldername = 'analyzed_data/gp1_calve/';
% dict.calving_gp1.save_fileprefix = 'ht_gp1_calve_';
% dict.calving_gp2.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_2.mat';
% dict.calving_gp2.save_foldername = 'analyzed_data/gp2_calve/';
% dict.calving_gp2.save_fileprefix = 'ht_gp2_calve_';
% dict.calving_gp3.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_3.mat';
% dict.calving_gp3.save_foldername = 'analyzed_data/gp3_calve/';
% dict.calving_gp3.save_fileprefix = 'ht_gp3_calve_';
% dict.calving_gp4.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_4.mat';
% dict.calving_gp4.save_foldername = 'analyzed_data/gp4_calve/';
% dict.calving_gp4.save_fileprefix = 'ht_gp4_calve_';
% 
% dict.calving_mu_gp1.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_1.mat';
% dict.calving_mu_gp1.save_foldername = 'analyzed_data/gp1_mu_calve';
% dict.calving_mu_gp1.save_fileprefix = 'ht_gp1_mu_calve_';
% dict.calving_mu_gp2.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_2.mat';
% dict.calving_mu_gp2.save_foldername = 'analyzed_data/gp2_mu_calve';
% dict.calving_mu_gp2.save_fileprefix = 'ht_gp2_mu_calve_';
% dict.calving_mu_gp3.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_3.mat';
% dict.calving_mu_gp3.save_foldername = 'analyzed_data/gp3_mu_calve';
% dict.calving_mu_gp3.save_fileprefix = 'ht_gp3_mu_calve_';
% dict.calving_mu_gp4.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_4.mat';
% dict.calving_mu_gp4.save_foldername = 'analyzed_data/gp4_mu_calve';
% dict.calving_mu_gp4.save_fileprefix = 'ht_gp4_mu_calve_';
%% save data from model outputs into .mat
% here we are saving 1. the sampled control point distance from the inflow
% boundary, 2. the time axis, and 3. the thickness at the sampled control
% points
% specify the experiment that you want to extract data from
mddict = dict.calving_pulse8;

modelname_calving = mddict.modelname;
save_foldername  = mddict.save_foldername;
save_fileprefix = mddict.save_fileprefix;
ds_i = sample_interval/ds;

% sorting by natsortfiles
folder_dir = natsortfiles(dir([pwd '/' foldername]));

plot_idx = 0;
md_count = 0;
for i = 1:size(folder_dir,1)
    % skip the irrelevant ones
    if ~strcmp(folder_dir(i).name(1), 'm')
        continue
    else
        md_count = md_count + 1;
        fullname = [folder_dir(i).name '/' modelname_calving];
        % load the ISSM model
        calve_md = load([folder_dir(i).folder '/' fullname]).md;
        nt = size(calve_md.results.TransientSolution,2);

        % Get the grounding line and calving front location sequence
        gl_locs = zeros(1,nt);
        front_locs = zeros(1,nt);
        for k = 1:nt
            front_mask = calve_md.results.TransientSolution(k).MaskIceLevelset;
            gl_mask = calve_md.results.TransientSolution(k).MaskOceanLevelset;
            front_locs(k) = locate_calvingfront(calve_md, front_mask);
            gl_locs(k) = locate_groundingline(calve_md, gl_mask);
        end
        
        % Get the thickness, control point location, and time axis
        % geometry parameters
        Lx = max(calve_md.mesh.x);
        Ly = max(calve_md.mesh.y);
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        % save the mesh elements and (x,y)
        mesh_elements = calve_md.mesh.elements;
        mesh_x = calve_md.mesh.x;
        mesh_y = calve_md.mesh.y;
        % retrieve the last position of the ice front; spclevelset ~= 0
        pos = find(calve_md.levelset.spclevelset(1:end-1,end) < 1 & calve_md.levelset.spclevelset(1:end-1,end) > -1);
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
        calve_results = struct2table(calve_md.results.TransientSolution);
        empty_md = calve_md;
        empty_md.results.TransientSolution = [];
        clear calve_md

        % initialize space to store h(t)
        thalweg_sample_ht = [];

        % get all h(t) data
        for j = 1:nt
            calve_surface = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
                calve_results.Thickness{j},...
                x, y, NaN);
            surface_profile  = calve_surface(mid_i,sample_i);
            thalweg_sample_ht = [thalweg_sample_ht; surface_profile];
        end
        thalweg_sample_ht = thalweg_sample_ht - thalweg_sample_ht(1,:);
        % crop the one extra time step in calving simulation
        thalweg_sample_ht = thalweg_sample_ht(1:end-1,:);

        % time vector
        time = calve_results.time;
        % shift time to start at 0
        time = time(1:end-1);
        time = time - time(1);

        % save data
        ht_data.h = thalweg_sample_ht;
        ht_data.t = time;
        ht_data.x = sample_x;
        ht_data.gl = gl_locs;
        ht_data.front = front_locs;
        filename = [save_foldername, save_fileprefix, folder_dir(i).name,'.mat'];
        save(filename,'ht_data')

        % print
        disp(['Model ',fullname,' is saved!'])
    end
end