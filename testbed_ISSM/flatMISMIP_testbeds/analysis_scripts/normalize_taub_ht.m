%% find if d(deltaH)/d(delta_Tau_b) is similar across all glaciers
% parameters
pulse_type = 'Pulse'; % types: "Diffu","Pulse"
geom_type = 'deep'; % types: "deep", "shallow"
expt_type = 'no_mu'; % types: "no_mu" (no mass unloading) and "mu" (with mass unloading)

% constants
ds = 50;
gauss_xloc = 3.2e4;
gauss_width_ratio = 0.08;
wid_eff_factor = 4; % 4 times the sigma (in the gaussian width) 

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
switch expt_type
    case 'no_mu'
        expt_name = ['MISMIP_yangTransient_Calving_',pulse_type,'GaussianPerturb_8.mat'];
        ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
    case 'mu'
        expt_name = ['MISMIP_yangTransient_Calving_MassUnloading_',pulse_type,'GaussianPerturb_8.mat'];
        ctrl_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
    otherwise
        error('Experiment type unrecognized!')
end
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

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

n_simu = size(folder_dir_groups{geom_i}, 1);
slopes = zeros(1,n_simu);

% iterate over deep or shallow GL models
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    results_tbl_expt = struct2table(md_expt.results.TransientSolution);
    results_tbl_ctrl = struct2table(md_ctrl.results.TransientSolution);
    % get model information
    modelname = md_expt.miscellaneous.name;
    [W, GL, FC] = parse_modelname(modelname);
    % isolate the delta H from localized basal perturbation
    deltaH = [results_tbl_expt.Thickness{:}] - [results_tbl_ctrl.Thickness{:}];
    deltaMask = [results_tbl_expt.MaskOceanLevelset{:}] - [results_tbl_ctrl.MaskOceanLevelset{:}];
    deltaH_cell = num2cell(deltaH,1);
    [deltaH_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
    [mask_grid, ~, ~] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, results_tbl_ctrl.MaskIceLevelset, ds);
    % permute
    deltaH_grid = permute(deltaH_grid,[2,3,1]);
    
    % apply the last mask to data at all timesteps
    mask = mask_grid(end,:,:);
    for i = 1:size(deltaH_grid,3)
        md_temp = deltaH_grid(:,:,i);
        md_temp(mask >=0) = 0;
        deltaH_grid(:,:,i) = md_temp;
    end
    % similarly, find delta_tau_b from expt-ctrl, map to grid, and mask
    % calculate basal drag
    vel = [results_tbl_expt.Vel{:}]/md_expt.constants.yts;
    % next step is necessary if the experiment is the combined experiment
    [lia, lib] = ismembertol(md_expt.friction.C(end,:),results_tbl_expt.time, 1e-5);
    idx_bigtimestep = sort(find(lib~=0));
    if length(idx_bigtimestep)~=length(results_tbl_expt.time)
        error('Corresponding shorter time axis is not found propoerly!')
    end
    C_cut = md_expt.friction.C(1:end-1,idx_bigtimestep);
    tau_b = C_cut.^2.*vel;
    % get vel and tau_b from ctrl
    vel_ctrl = [results_tbl_ctrl.Vel{:}]/md_ctrl.constants.yts;
    if size(md_ctrl.friction.C,2) == 1 % the "control" has no mass unloading
        C_mtx = repmat(md_ctrl.friction.C, 1,size(tau_b,2));
        vel_ctrl = [results_tbl_ctrl.Vel{:}]/md_ctrl.constants.yts;
        tau_b_ctrl = C_mtx.^2.*vel_ctrl;

    else % the "control" has mass unloading
        % similar to before, we search for the corresponding dt = 0.1 yr
        % timestep and remove the rest
        [lia, lib] = ismembertol(md_ctrl.friction.C(end,:),results_tbl_expt.time, 1e-5);
        idx_bigtimestep = sort(find(lib~=0));
        if length(idx_bigtimestep)~=length(results_tbl_expt.time)
            error('Corresponding shorter time axis is not found propoerly!')
        end
        C_cut_ctrl = md_ctrl.friction.C(1:end-1,idx_bigtimestep);
        tau_b_ctrl = C_cut_ctrl.^2.*vel_ctrl;
    end
        
    % get change in basal drag
    delta_taub = tau_b - tau_b_ctrl;
    delta_taub_cell = num2cell(delta_taub,1);
    % convert to grid
    [delta_taub_grid, x, y] = mesh_to_grid_overtime(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, delta_taub_cell, ds);
    % permute
    delta_taub_grid = permute(delta_taub_grid,[2,3,1]);
    
    % apply the last mask to data at all timesteps
    mask = mask_grid(end,:,:);
    for i = 1:size(delta_taub_grid,3)
        md_temp = delta_taub_grid(:,:,i);
        md_temp(mask >=0) = 0;
        delta_taub_grid(:,:,i) = md_temp;
    end

    % we truncate the glaciers (10 km up and downstream of the
    % perturbation); we also limit the search to a 5 km band (2.5km up and
    % 2.5 km down) around the center line
    loc_up = gauss_xloc - 10000;
    loc_down = gauss_xloc + 10000;
    loc_left = max(y)/2 - 1000;
    loc_right = max(y)/2 + 1000;
    
    x_i = find(x <= loc_down & x >= loc_up); 
    y_i = find(y <= loc_right & y >= loc_left);
    x_crop = x(x_i);
    y_crop = y(y_i);
    delta_taub_gcrop = delta_taub_grid(y_i,x_i,:);
    deltaH_gcrop = deltaH_grid(y_i,x_i,:);
    % find the loc of maximum and minimum amplitude in H(t)
    [~,max_i] = max(deltaH_gcrop,[],'all');
    [~,min_i] = min(deltaH_gcrop,[],'all');
    [max_yi, max_xi,~] = ind2sub(size(deltaH_gcrop),max_i);
    [min_yi, min_xi,min_ti] = ind2sub(size(deltaH_gcrop),min_i);
    % find the loc of most decreased basa drag
    [~,min_k] = min(delta_taub_gcrop,[],'all');
    [min_yk, min_xk, min_tk] = ind2sub(size(delta_taub_gcrop),min_k);

    % extract the timeseries from those points
    Ht_min = squeeze(deltaH_gcrop(min_yi,min_xi,:)); 
    Ht_max = squeeze(deltaH_gcrop(max_yi,max_xi,:));
    taub_min = squeeze(delta_taub_gcrop(min_yk,min_xk,:));

    % draw on the graph
    figure('Position',[100,100,900,500]); 
    tiledlayout(4,1,'TileSpacing','none')
    nexttile % H(t)
    imagesc(x_crop,y_crop,deltaH_gcrop(:,:,min_ti)); hold on; clim([-10,10]); colorbar
    scatter(x_crop(max_xi),y_crop(max_yi),40,'filled','red'); hold on
    scatter(x_crop(min_xi),y_crop(min_yi),40,'filled','yellow'); hold off
    title('delta H')
    nexttile % basa drag
    imagesc(x_crop,y_crop,delta_taub_gcrop(:,:,min_tk)); hold on; clim([-3e4, 3e4]); colorbar
    scatter(x_crop(min_xk),y_crop(min_yk),40,'filled','red'); 
    title('basal drag')
    nexttile([2,1])
    yyaxis left
    plot(Ht_min);hold on;
    plot(Ht_max);hold on; 
    ylabel('elevation change (m)')
    yyaxis right; 
    plot(taub_min)
    ylabel('basal drag change (Pa)')
    legend(["min delta H","max delta H","min delta tau_b"])
    % save to files for reference
    filename = ['plots/locate_peak/',modelname(9:end),'_',pulse_type,'_',geom_type,'_',expt_type,'.pdf'];
    exportgraphics(gcf,filename,'ContentType','vector')

    % find peaks 
    [pks_ht, locs_ht] = findpeaks(-1*Ht_min(50:end));
    [pks_taub, locs_taub] = findpeaks(-1*taub_min(50:end));
    % find slope
    p = polyfit(pks_taub, pks_ht,1);
    slopes(j) = p(1);

    %% detrend, locate the peaks, find mean peak amplitude
    Ht_min_c = Ht_min(50:end);
    taub_min_c = taub_min(50:end);
    figure;
    plot(Ht_min);hold on
    plot(detrend(Ht_min,15))
    legend(["original","detrended"])
    
    n = 15;
    [pks_ht, locs_ht] = findpeaks(detrend(Ht_min_c,n));
    [pks_ht_neg, ~] = findpeaks(-1*detrend(Ht_min_c,n));
    pks_ht_neg = -1*sort(pks_ht_neg,'descend');
    pks_ht = sort(pks_ht,'descend');
    % keep the largest 8 (8 pulses in total)
    pks_ht = pks_ht(1:8);
    pks_ht_neg = pks_ht_neg(1:8);
    mean_amp = mean(pks_ht) - mean(pks_ht_neg);
    

    disp(['Model ',modelname(9:end),' is completed!'])
end