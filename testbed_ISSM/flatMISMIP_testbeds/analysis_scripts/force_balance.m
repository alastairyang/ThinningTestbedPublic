%% force balance
% We find the temporal change in force balance structure
%% Main script
geom_type = "deep"; % options: "deep" or "shallow"
ds = 400; % grid size for the regular grid
sampled_ti = 1:10:240; % sampled time index

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
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
driving_S_all = cell(2, 260, n_simu);
longi_grad_all = cell(2, 260, n_simu);
later_grad_all = cell(2, 260, n_simu);
basal_R_all = cell(2, 260, n_simu);
gl_x_all = cell(2, 260, n_simu);
front_x_all = cell(2, 260, n_simu);
Ws_md  = zeros(1,n_simu);
GLs_md = zeros(1,n_simu);
FCs_md = zeros(1,n_simu);
dH_max_expt = zeros(1,n_simu);
dH_max_ctrl = zeros(1,n_simu);
dH_sum_expt = zeros(1,n_simu);
dH_sum_ctrl = zeros(1,n_simu);
gl_expt = zeros(1,n_simu);
gl_ctrl = zeros(1,n_simu);

% get both the force balance components and 
for j = 1:n_simu
    % read the model
    group = folder_dir_groups{geom_i};
    md_expt = load([group.folder{j},'/', group.name{j}, '/', expt_name]).md;
    md_ctrl = load([group.folder{j},'/', group.name{j}, '/', ctrl_name]).md;
    index = md_expt.mesh.elements;
    % get the maximum thinning along the centerline
    % experiment
    final_icemask = md_expt.results.TransientSolution(end).MaskIceLevelset;
    final_H = md_expt.results.TransientSolution(end).Thickness; final_H(final_icemask>0) = 0;
    first_H = md_expt.results.TransientSolution(1).Thickness;   first_H(final_icemask>0) = 0;
    deltaH = mesh_to_grid(md_expt.mesh.elements, md_expt.mesh.x, md_expt.mesh.y, first_H-final_H, ds);
    yi_mid = floor(size(deltaH,1)/2);
    dH_sum_expt(j) = sum(deltaH(yi_mid,:))*ds;
    dH_max_expt(j) = max(deltaH(yi_mid,:),[],'all');
    % control
    final_icemask = md_ctrl.results.TransientSolution(end).MaskIceLevelset;
    final_H = md_ctrl.results.TransientSolution(end).Thickness; final_H(final_icemask>0) = 0;
    first_H = md_ctrl.results.TransientSolution(1).Thickness;   first_H(final_icemask>0) = 0;
    deltaH = mesh_to_grid(md_ctrl.mesh.elements, md_ctrl.mesh.x, md_ctrl.mesh.y, first_H-final_H, ds);
    yi_mid = floor(size(deltaH,1)/2);
    dH_sum_ctrl(j) = sum(deltaH(yi_mid,:))*ds;
    dH_max_ctrl(j) = max(deltaH(yi_mid,:),[],'all');
    % record model info
    [Ws_md(j), GLs_md(j), FCs_md(j)] = parse_modelname(md_expt.miscellaneous.name);

    % get the maximum grounding line retreat dist
    % first experiment
    final_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_expt,md_expt.results.TransientSolution(1).MaskOceanLevelset);
    gl_expt(j) = abs(final_gl - first_gl);
    % then control
    final_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(end).MaskOceanLevelset);
    first_gl = locate_groundingline(md_ctrl,md_ctrl.results.TransientSolution(1).MaskOceanLevelset);
    gl_ctrl(j) = abs(final_gl - first_gl);

    % get force balance components
    for ti = sampled_ti
        smooth_L = 2000;
        [driving_S_expt, basal_R_expt, longi_grad_expt, later_grad_expt, x, y, gl_x_expt, front_x_expt] = calc_force_balance(md_expt,ti,smooth_L);
        [driving_S_ctrl, basal_R_ctrl, longi_grad_ctrl, later_grad_ctrl, ~, ~, gl_x_ctrl, front_x_ctrl] = calc_force_balance(md_ctrl,ti,smooth_L);
        % save
        % experiment
        driving_S_all{1,ti,j} = driving_S_expt; 
        longi_grad_all{1,ti,j} = longi_grad_expt;
        later_grad_all{1,ti,j} = later_grad_expt;
        basal_R_all{1,ti,j} = basal_R_expt;
        gl_x_all{1,ti,j} = gl_x_expt;
        front_x_all{1,ti,j} = front_x_expt;
        % control
        driving_S_all{2,ti,j} = driving_S_ctrl; 
        longi_grad_all{2,ti,j} = longi_grad_ctrl;
        later_grad_all{2,ti,j} = later_grad_ctrl;
        basal_R_all{2,ti,j} = basal_R_ctrl;
        gl_x_all{2,ti,j} = gl_x_ctrl;
        front_x_all{2,ti,j} = front_x_ctrl;

    end
    disp(['model ',num2str(j), ' is completed!'])
end

%% Get total stress balance
runme_params = readtable('runme_param.csv');
front_x = runme_params.terminus0_x;
total_Rs = zeros(2,n_simu);

for ri = 1:2 % first expt, then control
    for j = 1:n_simu % simulations
        ti = sampled_ti(1);
        % get a referential resistive stress estimate at the first timestep
        yi_mid = floor(size(basal_R_all{ri,ti,j},1)/2);
        sample_wid = floor(Ws_md(j)*0.6/ds/2);
        yi_central = yi_mid-sample_wid:yi_mid+sample_wid;
        % calculate the resistive stress
        later_grad_central = later_grad_all{ri,ti,j}(yi_central,:);
        later_grad_mid = mean(later_grad_central,1);
        basal_R_central = basal_R_all{ri,ti,j}(yi_mid,:);
        basal_R_mid = mean(basal_R_central,1);

        % first: we find the integrated resistive stress
        % from gl(end time) to the calving front when basal
        % shear stress is still in place
        %gl_x = gl_x_all{ri,sampled_ti(end),j};
        gl_x = 0;
        front_x = front_x_all{ri,sampled_ti(end),j};
        sample_xi = find(gl_x < x & x < front_x);
        % integrate with trapezoid method
        later_int = trapz(later_grad_mid(sample_xi))*ds;
        basal_int = trapz(basal_R_mid(sample_xi))*ds;
        total_R_ref = later_int+basal_int;

        % at the last timestep, get an estimate of loss of total resistive
        % stress
        ti = sampled_ti(end);
        %gl_x = gl_x_all{ri,ti,j};
        gl_x = 0;
        front_x = front_x_all{ri,ti,j};
        sample_xi = find(gl_x < x & x < front_x);
        later_grad_sample = later_grad_mid(sample_xi);
        later_int = trapz(later_grad_sample)*ds;
        total_R_last = later_int;

        % find the difference
        total_Rs(ri,j) = total_R_ref - total_R_last;
        clear total_R_ref total_R_last

    end
               
end
%exportgraphics(gcf,'plots/taub_taud_frac.png','Resolution',300)

%% Plotting maximum thinning and total resistive stress loss 
Ws_symb = [40,100,260];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;
figure('Position',[100,100,1100,600]);
tiledlayout(1,3,'TileSpacing','compact')

nexttile % total thinning vs GL retreat
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(gl_expt(j)/1e3, dH_sum_expt(j),W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(gl_ctrl(j)/1e3, dH_sum_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
ylabel('Total thinning (m^2)')
xlabel('Grounding line retreat (km)')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Aria';

nexttile % Max thinning vs GL retreat
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(gl_expt(j)/1e3, dH_max_expt(j), W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(gl_ctrl(j)/1e3, dH_max_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
ylabel('Maximum thinning (m)')
xlabel('Grounding line retreat (km)')
set(gca,'XTick',[2,6,10,14]);
ax = gca; ax.FontSize = 14; ax.FontName = 'Aria';

% nexttile
% for j = 1:n_simu
%     W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
%     GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
%     FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
%     % plot the experiment
%     scatter(dH_sum_expt(j),total_Rs(1,j)/1e6, W_symb,FC_symb,'filled',GL_symb)
%     hold on
%     % plot the control
%     scatter(dH_sum_ctrl(j), total_Rs(2,j)/1e6, W_symb,FC_symb,GL_symb);
%     hold on
% end
% xlabel('Total thinning (m^2)','FontName','Aria','FontSize',14)
% ylabel('Total resistive stress loss (MPa)','FontName','Aria','FontSize',14)

nexttile % Total stress loss vs max thinning
for j = 1:n_simu
    W_symb = Ws_symb(Ws_md(j)==Ws); % marker size
    GL_symb = GLs_symb(GLs_md(j)==GLs); % marker type (square is shallow; circle is deep)
    FC_symb = FCs_symb(FCs_md(j)==FCs,:); % color
    % plot the experiment
    scatter(total_Rs(1,j)/1e9, dH_max_expt(j), W_symb,FC_symb,'filled',GL_symb)
    hold on
    % plot the control
    scatter(total_Rs(2,j)/1e9, dH_max_ctrl(j), W_symb,FC_symb,GL_symb);
    hold on
end
xlim([0,100])
ylabel('Maximum thinning (m)','FontName','Aria','FontSize',14)
xlabel('Total resistive stress loss (GPa)','FontName','Aria','FontSize',14)
ax = gca; ax.FontSize = 14; ax.FontName = 'Aria';

exportgraphics(gcf,'plots/two_thinning_proxies.png','Resolution',600)


%% plot evolution of the absolute force component
figure('Position',[100,100,1200,600]);
tiledlayout(3,3,'TileSpacing','none')
for j = 1:n_simu
    nexttile
    colororder(cool(length(sampled_ti)))
    for ti = sampled_ti 
        % plot longitudinal component
        yi_mid = floor(size(basal_R_all{ti,j},1)/2);
        basal_R_mid   = basal_R_all{ti,j}(yi_mid,:);
        plot(x/1000, basal_R_mid)
        hold on;
        ylim([0,3.1e5])
    end

    % modify the tiled plot appearance
    if j == 4
        ylabel('$\tau_b (Pa)$','Interpreter','latex','FontSize',20)    
    end
    if j == 8
        xlabel('Along flow distance (km)','Interpreter','latex','FontSize',20)
    end
    if ismember(j, [2,3,5,6,8,9])
        set(gca,'ytick',[]);
    else
        set(gca,'ytick',[1e5,2e5, 3e5])
    end
    if ismember(j, [1,2,3,4,5,6])
        set(gca,'xtick',[]); 
    else
        set(gca,'xtick',[10,20,30,40])
    end

end
exportgraphics(gcf,'plots/taub.png','Resolution',300)
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

