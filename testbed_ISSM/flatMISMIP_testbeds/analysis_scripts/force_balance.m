%% force balance
% here we find the change in fractional basal resistive stress at the end
% of the perturbation. 
% Unlike th old code, we right now ignore the longitudinal and the lateral
% resistive stress; also we look at all the models rather than 2
geom_type = "deep"; % options: "deep" or "shallow"
ds = 50; % grid size for the regular grid
front_x = 56650; % terminus distance to x = 0

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
md_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
% get all model foldernames
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);

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
n_simu = size(folder_dir_groups{shallowGL_i}, 1);

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown depth specification!')
end

% save the force balance field data
driving_S_all = cell(2, n_simu);
longi_grad_all = cell(2, n_simu);
later_grad_all = cell(2, n_simu);
for j = 9
    % read the model
    group = folder_dir_groups{geom_i};
    md = load([group.folder{j},'/', group.name{j}, '/', md_name]).md;
    index = md.mesh.elements;
    %compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
    [alpha, beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);
    summation=[1;1;1];
    
    nt = size(md.results.TransientSolution,2);
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 50;
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);

    % get H from vertices to elements
    % timesteps we look at: first and last
    timesteps = [1, size(md.results.TransientSolution,2)];
    for ti = 1:length(timesteps)
        t = timesteps(ti);
        % ice thickness
        H = md.results.TransientSolution(t).Thickness;
        H_list = H(index);
        H_list = mean(H_list,2);

        % basal stress; get onto elements
        if size(md.friction.C,2) == 1
            % no sliding law coefficient change
            bs = md.friction.C.^2.*md.results.TransientSolution(t).Vel/md.constants.yts;
        else
            % mass unloading experiment
            bs = md.friction.C(1:end-1,t).^2.*md.results.TransientSolution(t).Vel/md.constants.yts;
        end
        bs_list = bs(index);
        basal_R = mean(bs_list,2);
        % lateral
        n = 3;
        Rxx = md.materials.rheology_B.*md.results.TransientSolution(t).StrainRateeffective.^(1/n-1).*(2*md.results.TransientSolution(t).StrainRatexx +   md.results.TransientSolution(t).StrainRateyy);
        Ryy = md.materials.rheology_B.*md.results.TransientSolution(t).StrainRateeffective.^(1/n-1).*(  md.results.TransientSolution(t).StrainRatexx + 2*md.results.TransientSolution(t).StrainRateyy);
        Rxy = md.materials.rheology_B.*md.results.TransientSolution(t).StrainRateeffective.^(1/n-1).*md.results.TransientSolution(t).StrainRatexy;
        
        % driving stress
        driving_S = drivingstress_from_results(md, t);
        
        % make floating part NaN
        mask = md.results.TransientSolution(t).MaskOceanLevelset;
        %basal_R(mask<0) = nan;

        % make grid
        Lx = max(md.mesh.x);
        Ly = max(md.mesh.y);
        ds = 50;
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        % interp onto grids
        [driving_S, ~, ~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, driving_S, ds);
        [basal_R, ~, ~]   = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, basal_R, ds);
        % smooth (along-flow direction) the driving stress at length scale
        % of 1 km
        n_ds = 1000/ds;
        driving_S_smooth = movmean(driving_S, n_ds, 2); % smooth() averages along each column vector
        %driving_S_smooth = imgaussfilt(driving_S, n_ds);

        % calculate longitudinal and lateral stress gradient
        % we notice that taking derivative with the raw Rxx and Rxy gives
        % very noisy results. To mitigate, we first interpolate onto a
        % coarse grid and then interpolate to a finer grid
        % with cubic spline, then finally take derivative with five-point stencil.
        ds = 200;
        [Rxxgrid,~,~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, Rxx, ds);
        [Rxygrid,~,~] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, Rxy, ds);
        [Hgrid,x,y] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, H, ds);
        ds = 50; % back to higher resolution
        % get the new meshgrid
        [~,xq,yq] = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, H, ds);
        [X,Y] = meshgrid(x,y);
        [Xq, Yq] = meshgrid(xq, yq);
        Rxxgrid_interp = interp2(X,Y,Rxxgrid,Xq,Yq,'spline');
        Rxygrid_interp = interp2(X,Y,Rxygrid,Xq,Yq,'spline');
        Hgrid_interp   = interp2(X,Y,Hgrid,Xq,Yq, 'spine');
        % take spatial derivative with five-point stencile for higher
        % accuracy
        longi_grad = -1*five_point_stencil(Rxxgrid_interp.*Hgrid_interp, ds, 2);
        later_grad = -1*five_point_stencil(Rxygrid_interp.*Hgrid_interp, ds, 1);
        
        % save
        % first row: first time step; second row: last time step
        driving_S_all{ti,j} = basal_R./driving_S_smooth; 
        longi_grad_all{ti,j} = longi_grad;
        later_grad_all{ti,j} = later_grad;
    end
end

%% plot t = 1 and t = end
% t = 1
figure('Position',[100,100,800,300]);
tiledlayout(3,3,'TileSpacing','none')
for i = 1:n_simu
    nexttile
    imagesc(fb_ratio{1,i});
    clim([0,1])
end

figure('Position',[600,100,800,300]);
tiledlayout(3,3,'TileSpacing','none')
for i = 1:n_simu
    nexttile
    imagesc(fb_ratio{2,i});
    clim([0,1])
end
%% APPENDIX: Old force balance script
md1_name = "model_W5000_GL0_FC120000";
md2_name = "model_W5000_GL0_FC30000";
model_type = "MISMIP_yangTransient_Calving_MassUnloading.mat";
md1 = load("long_models_yang/" + md1_name + "/" + model_type).md;
md2 = load("long_models_yang/" + md2_name + "/" + model_type).md;
mds = [md1, md2];

figure;
for md_i = 1:length(mds)
    md = mds(md_i);
    index = md.mesh.elements;
    %compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
    [alpha, beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);
    summation=[1;1;1];
    
    nt = size(md.results.TransientSolution,2);
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 50;
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);
    
    color_length = nt;
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
        linspace(red(2),sth(2),color_length)',...
        linspace(red(3),sth(3),color_length)'];
    fb_ratio_last = 0;
    iter_count = 0;

    % iterate over time
    for i = 30:5:nt
        iter_count = iter_count + 1;
        tauxx = md.results.TransientSolution(i).DeviatoricStressxx;
        tauxy = md.results.TransientSolution(i).DeviatoricStressxy;
        tauxxlist=tauxx(index);
        tauxylist=tauxy(index);
        % get H from vertices to elements
        H = md.results.TransientSolution(i).Thickness;
        H_list = H(index);
        H_list = mean(H_list,2);
        % find directional derivative along x, y
        dtauxxdx=(tauxxlist.*H_list.*alpha)*summation;
        dtauxxdy=(tauxxlist.*H_list.*beta)*summation;
        dtauxydx=(tauxylist.*H_list.*alpha)*summation;
        dtauxydy=(tauxylist.*H_list.*beta)*summation;
        % basal stress; get onto elements
        if size(md.friction.C,2) == 1
            % no sliding law coefficient change
            bs = md.friction.C.^2.*md.results.TransientSolution(i).Vel/md.constants.yts;
        else
            % mass unloading experiment
            bs = md.friction.C(1:end-1,i).^2.*md.results.TransientSolution(i).Vel/md.constants.yts;
        end
        bs_list = bs(index);
        bs = mean(bs_list,2);
        % driving stress
        ds = drivingstress_from_results(md, i);
        
        time = md.results.TransientSolution(i).time;
        plot_title = [md.miscellaneous.name, ', time = ', num2str(time)];
        mask = md.results.TransientSolution(i).MaskOceanLevelset;
        mask = mean(mask(index),2);
    
        % force balance: longitudinal-lateral stress gradient / driving
        % stress
        % The part of d(tauxx)/dx that contributes to the driving stress,
        % we single it out, remove from resistive calculation, and add to
        % the driving stress
        ds_dtauxxdx = zeros(size(dtauxxdx));
        ds_dtauxxdx(dtauxxdx>0) = dtauxxdx(dtauxxdx>0);
        dtauxxdx(dtauxxdx>0) = 0;
        resistive = -(dtauxxdx + dtauxydy);
        fb_ratio = resistive./(ds + ds_dtauxxdx);
        
        fb_ratio(mask<0) = nan;
        
%        subplot(1,length(mds),md_i)

        Lx = max(md.mesh.x);
        Ly = max(md.mesh.y);
        ds = 250;
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        
        mask = md.results.TransientSolution(i).MaskOceanLevelset;
        field_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            fb_ratio, x, y, NaN);
        field_profile = field_grid(mid_i,:);
        field_prof_smooth = smooth(field_profile(~isnan(field_profile)),15);
        field_profile(~isnan(field_profile)) = field_prof_smooth;
%         delta_field_profile = field_profile - field_profile_last;
%         % update field_profile_last
%         field_profile_last = field_profile;
%         % look for only 15 km behind the grounding line
        if iter_count == 1
            continue 
        else

%             gl_x = locate_groundingline(md,mask);
%             x_keep = find(thalweg_x > gl_x-1.5e4 & thalweg_x < gl_x);
%             % truncate
%             field_profile = field_profile(x_keep);
%             thalweg_x = thalweg_x(x_keep);
%             thalweg_x = thalweg_x - max(thalweg_x);
%         % plot
%             plot(thalweg_x, field_profile, Color=colors_p(i,:));hold on;
%             title(plot_title)
%             ylim([0,1])
%             pause(0.1)
            gif
            plotmodel(md,'data',fb_ratio - fb_ratio_last,'caxis',[0,0.2],'mask',mask)
            fb_ratio_last = fb_ratio;
            pause(0.2)
        end
        %plotmodel(md,'data',-dtauxxdx,'caxis',[0,2e4])
    end
end

