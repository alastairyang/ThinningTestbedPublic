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

dict.calving_smw.modelname = 'MISMIP_yangTransient_Calving_SMweakening.mat';
dict.calving_smw.save_foldername = 'analyzed_data/smw_calve/';
dict.calving_smw.save_fileprefix = 'ht_smw_calve_';

dict.calving_mu.modelname = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
dict.calving_mu.save_foldername = 'analyzed_data/mu_calve/';
dict.calving_mu.save_fileprefix = 'ht_mu_calve_';

dict.calving_gp1.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_1.mat';
dict.calving_gp1.save_foldername = 'analyzed_data/gp1_calve/';
dict.calving_gp1.save_fileprefix = 'ht_gp1_calve_';
dict.calving_gp2.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_2.mat';
dict.calving_gp2.save_foldername = 'analyzed_data/gp2_calve/';
dict.calving_gp2.save_fileprefix = 'ht_gp2_calve_';
dict.calving_gp3.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_3.mat';
dict.calving_gp3.save_foldername = 'analyzed_data/gp3_calve/';
dict.calving_gp3.save_fileprefix = 'ht_gp3_calve_';
dict.calving_gp4.modelname = 'MISMIP_yangTransient_Calving_GaussianPerturb_4.mat';
dict.calving_gp4.save_foldername = 'analyzed_data/gp4_calve/';
dict.calving_gp4.save_fileprefix = 'ht_gp4_calve_';

dict.calving_mu_gp1.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_1.mat';
dict.calving_mu_gp1.save_foldername = 'analyzed_data/gp1_mu_calve';
dict.calving_mu_gp1.save_fileprefix = 'ht_gp1_mu_calve_';
dict.calving_mu_gp2.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_2.mat';
dict.calving_mu_gp2.save_foldername = 'analyzed_data/gp2_mu_calve';
dict.calving_mu_gp2.save_fileprefix = 'ht_gp2_mu_calve_';
dict.calving_mu_gp3.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_3.mat';
dict.calving_mu_gp3.save_foldername = 'analyzed_data/gp3_mu_calve';
dict.calving_mu_gp3.save_fileprefix = 'ht_gp3_mu_calve_';
dict.calving_mu_gp4.modelname = 'MISMIP_yangTransient_Calving_MassUnloading_GaussianPerturb_4.mat';
dict.calving_mu_gp4.save_foldername = 'analyzed_data/gp4_mu_calve';
dict.calving_mu_gp4.save_fileprefix = 'ht_gp4_mu_calve_';
%% save __ data from model outputs into .mat
% specify the experiment that you want to extract data from
mddict = dict.calving_gp4;

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
        calve_md = load([folder_dir(i).folder '/' fullname]).md;
        nt = size(calve_md.results.TransientSolution,2);

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
        filename = [save_foldername, save_fileprefix, folder_dir(i).name,'.mat'];
        save(filename,'ht_data')

        % print
        disp(['Model ',fullname,' is saved!'])
    end
end

%% time lag to reversal: calving only
reversal_year = 5+8; % 5 year stationary in the beginning and 8 accelerated retreat years.

foldername = 'analyzed_data/calve_only';
folder_prefix = 'ht_calve_';
folder_dir = natsortfiles(dir([pwd '/' foldername]));
folder_dir = struct2table(folder_dir);
% remove  '.' and '..'
bools = cellfun(@(s) ~strcmp(s(1),'.'), folder_dir.name);
folder_dir = folder_dir(bools,:);

plot_idx = 0;

ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

GLs = [0,400];
folder_dir_groups = cell(1,2);
% split the folder_dir into two groups, separated by grounding line depth
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(folder_dir,1),1);
    for j = 1:size(folder_dir.name)
        GL_bool(j) = compare_GLvalue(folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end
        

figure('Position',[100,100,700,700]);
for j = 1:length(GLs) % iterate over grounding line depths
    folder_dir = folder_dir_groups{j};
    md_count = 0;

    for i = 1:size(folder_dir,1)
        % skip the irrelevant ones
        md_count = md_count + 1;

        md = load(string(folder_dir.folder(i))+"/"+ string(folder_dir.name(i))).ht_data;
        modelname = folder_dir.name(i);
        modelname = modelname{1}(length(folder_prefix)+1:end-4);
        % need to coarsen / smooth the data, or the derivative produces wiggles
        smooth_window = 20;
        % to smooth data, we have to iterate over each flow line
        % reshape the matrix into a vector and smooth will create
        % jumps/shocks
        h_smooth = zeros(size(md.h));
        for line_i = 1:size(md.h,2)
            h_smooth(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end
        dh = diff(h_smooth,1,1);
        dt = md.t(2) - md.t(1);
        dhdt = dh/dt;

        % find the lag year to calving front reversal
        [min_dhdt, min_dhdt_i] = min(dhdt,[], 1);
        min_dhdt_year = md.t(min_dhdt_i);
        % create along flowline sampling distance
        distances = (1:sample_number)*sample_interval;

        subplot(3,3, md_count)
        % plot the year lag over distance to ice front
        % we can choose to skip the first few sample points by specifying%
        % start_i bigger than 1
        start_i = 3;
        distances = distances(start_i:end);
        min_dhdt_year = min_dhdt_year(start_i:end);
        min_dhdt = min_dhdt(start_i:end);
        scatter(distances/1000, min_dhdt_year - reversal_year,abs(min_dhdt)*5)
        if j == length(GLs)
            legend(["GL = "+string(num2str(GLs(1))), "GL = "+string(num2str(GLs(2)))]);
        end
        hold on
        ylim([0,13])
%         subplot_title = strrep(modelname(1:end-4), '_',', ');
%         title(subplot_title)

        if ~ismember(md_count, ylabel_i)
            set(gca, 'ytick', [])
        end
        if ~ismember(md_count, xlabel_i)
            set(gca, 'xtick', [])
        end

    end
    subplot(3,3,8)
    xlabel('Distance to ice front (km)','Interpreter','latex',FontSize=13)
    subplot(3,3,4)
    ylabel('Lag time (yr)','Interpreter','latex',FontSize=13)
end


% save to folder
saveas(gcf, 'plots/calve_timelag.pdf')


%% time lag to reversal: calving + mass unloading
reversal_year = 5+8; % 5 year stationary in the beginning and 8 accelerated retreat years.

foldername = 'analyzed_data/mu_calve';
folder_prefix = 'ht_mu_calve_';
folder_dir = natsortfiles(dir([pwd '/' foldername]));
folder_dir = struct2table(folder_dir);
% remove  '.' and '..'
bools = cellfun(@(s) ~strcmp(s(1),'.'), folder_dir.name);
folder_dir = folder_dir(bools,:);

plot_idx = 0;
ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

GLs = [0,400];
folder_dir_groups = cell(1,2);
% split the folder_dir into two groups grounding line depths
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(folder_dir,1),1);
    for j = 1:size(folder_dir.name)
        GL_bool(j) = compare_GLvalue(folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end

% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
% specify the scatter plot data symbols
% GL: circle vs dot; FC: color, light to dark; W: size of the symbol
Ws_symb = [50,70,90];
GLs_symb = ["square","o"];
FCs_symb = [162,174,235;107,128,232;32,66,232]/255;

figure('Position',[100,100,500,500]);
for j = 1:length(GLs) % iterate over grounding line depths
    folder_dir = folder_dir_groups{j};
    md_count = 0;

    for i = 1:size(folder_dir,1)
        md_count = md_count + 1;

        md = load(string(folder_dir.folder(i))+"/"+ string(folder_dir.name(i))).ht_data;
        modelname = folder_dir.name(i);
        modelname = modelname{1}(length(folder_prefix)+1:end-4);
        % get the model parameter
        [W, GL, FC] = parse_modelname(modelname);

        % need to coarsen the data, or the derivative produces wiggles
        smooth_window = 20;
        % to smooth data, we have to iterate over each flow line
        % reshape the matrix into a vector and smooth will create
        % jumps/shocks
        h_smooth = zeros(size(md.h));
        for line_i = 1:size(md.h,2)
            h_smooth(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end
        dh = diff(h_smooth,1,1);
        dt = md.t(2) - md.t(1);
        dhdt = dh/dt;

        % find the lag year to calving front reversal
        [min_dhdt, min_dhdt_i] = min(dhdt,[], 1);
        min_dhdt_year = md.t(min_dhdt_i);
        % create along flowline sampling distance
        distances = (1:sample_number)*sample_interval;

        %subplot(3,3, md_count)
        % plot the year lag
        % we can choose to skip the first few sample points by specifying%
        % start_i bigger than 1
        start_i = 6;
        distances = distances(start_i:end);
        min_dhdt_year = min_dhdt_year(start_i:end);
        min_dhdt = min_dhdt(start_i:end);
        %plot(distances/1000, min_dhdt_year - reversal_year)%,abs(min_dhdt)*5)
        if j == length(GLs)
            legend(["GL = "+string(num2str(GLs(1))), "GL = "+string(num2str(GLs(2)))]);
        end
        hold on
        ylim([0,13])
         subplot_title = strrep(modelname(1:end-4), '_',', ');

%         if ~ismember(md_count, ylabel_i)
%             set(gca, 'ytick', [])
%         end
%         if ~ismember(md_count, xlabel_i)
%             set(gca, 'xtick', [])
%         end

        % more compact illustration: 1 scatter plot
        % get the linear fit to lag time - distance 
        P = polyfit(distances/1000, min_dhdt_year - reversal_year, 1);
        slope = P(1); 
        % min, mean, max lag time
        mean_tau = mean(min_dhdt_year - reversal_year);
        min_tau = abs(min(min_dhdt_year - reversal_year)-mean_tau);
        max_tau = abs(max(min_dhdt_year - reversal_year)-mean_tau);
        % get the corresponding symbols for this scatter plot
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        scatter(slope, mean_tau, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        errorbar(slope, mean_tau, min_tau, max_tau,'k')
        hold on

    end
    ylim([0,12])
%     subplot(3,3,8)
%     xlabel('Distance to ice front (km)','Interpreter','latex',FontSize=13)
%     subplot(3,3,4)
%     ylabel('Lag time (yr)','Interpreter','latex',FontSize=13)
end


% save to folder
saveas(gcf, 'plots/mu_calve_timelag.pdf')

%% time lag to reversal: seasonal calving + mass unloading
reversal_year = 5+8; % 5 year stationary in the beginning and 8 accelerated retreat years.

foldername = 'analyzed_data/mu_seasoncalve';
folder_prefix = 'ht_mu_seasoncalve_';
folder_dir = natsortfiles(dir([pwd '/' foldername]));
folder_dir = struct2table(folder_dir);
% remove  '.' and '..'
bools = cellfun(@(s) ~strcmp(s(1),'.'), folder_dir.name);
folder_dir = folder_dir(bools,:);

plot_idx = 0;

ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

GLs = [0,400];
folder_dir_groups = cell(1,2);
% split the folder_dir into two groups, separated by grounding line depth
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(folder_dir,1),1);
    for j = 1:size(folder_dir.name)
        GL_bool(j) = compare_GLvalue(folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    folder_dir_groups{i} = folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end
        

figure('Position',[100,100,700,700]);
for j = 1:length(GLs) % iterate over grounding line depths
    folder_dir = folder_dir_groups{j};
    md_count = 0;

    for i = 1:size(folder_dir,1)
        % skip the irrelevant ones
        md_count = md_count + 1;

        md = load(string(folder_dir.folder(i))+"/"+ string(folder_dir.name(i))).ht_data;
        modelname = folder_dir.name(i);
        modelname = modelname{1}(length(folder_prefix)+1:end-4);
        % need to coarsen / smooth the data, or the derivative produces wiggles
%         t_coarsen = md.t(1):0.6:md.t(end);
%         h_coarsen = interp1(md.t, md.h, t_coarsen);
        smooth_window = 20;
        % to smooth data, we have to iterate over each flow line
        % reshape the matrix into a vector and smooth will create
        % jumps/shocks
        h_smooth = zeros(size(md.h));
        for line_i = 1:size(md.h,2)
            h_smooth(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end
        dh = diff(h_smooth,1,1);
        dt = md.t(2) - md.t(1);
        dhdt = dh/dt;

        % find the lag year to calving front reversal
        [min_dhdt, min_dhdt_i] = min(dhdt,[], 1);
        min_dhdt_year = md.t(min_dhdt_i);
        % create along flowline sampling distance
        distances = (1:sample_number)*sample_interval;

        subplot(3,3, md_count)
        % plot the year lag
        % we can choose to skip the first few sample points by specifying%
        % start_i bigger than 1
        start_i = 3;
        distances = distances(start_i:end);
        min_dhdt_year = min_dhdt_year(start_i:end);
        min_dhdt = min_dhdt(start_i:end);
        plot(distances/1000, min_dhdt_year - reversal_year)%abs(min_dhdt)*5)
        if j == length(GLs)
            legend(["GL = "+string(num2str(GLs(1))), "GL = "+string(num2str(GLs(2)))]);
        end
        %plot(t_coarsen(1:end-1), dhdt)
        hold on
        ylim([0,13])
         subplot_title = strrep(modelname(1:end-4), '_',', ');
         %title(subplot_title)

        if ~ismember(md_count, ylabel_i)
            set(gca, 'ytick', [])
        end
        if ~ismember(md_count, xlabel_i)
            set(gca, 'xtick', [])
        end

    end
    subplot(3,3,8)
    xlabel('Distance to ice front (km)','Interpreter','latex',FontSize=13)
    subplot(3,3,4)
    ylabel('Lag time (yr)','Interpreter','latex',FontSize=13)
end


% save to folder
saveas(gcf, 'plots/mu_seasoncalve_timelag.pdf')

%% functions
function bool = compare_GLvalue(str, val)
    filename_split = split(str, '_');
    initials = string(cellfun(@(s) s(1:2), filename_split, 'UniformOutput', false));
    GLvalue = filename_split(strcmp('GL', initials));
    GLvalue = str2double(GLvalue{1}(3:end));

    bool = GLvalue == val;
end


