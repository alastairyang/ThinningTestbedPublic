%% Summarize average thinning rate and attenuation distance in one scatter grid plot
% effective pressure experiment

% Attenuation distance is defined as to distance to the ice front where the
% 63% of total culmulative thinning has been registered.

%% Main script
% load model parameter
sim_params = readtable('runme_param.csv');

% parameters
ds = 100; % meshgrid spacing
t_total = sim_params.perturb_duration; % total number of years

% input
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));
ctrl_name = 'MISMIP_yangTransient_CalvingOnly.mat';
expt_name = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
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

% pre-allocate
n_simu = size(folder_dir_groups{1}, 1);
dH = cell(2,2,n_simu);
maxdHdt = zeros(2,n_simu);
dL = zeros(2,n_simu);
Ws_a = zeros(2,9); FCs_a = zeros(2,9); GLs_a = zeros(2,9);

% acquire average thinning rate and attenuation distance
md_types = ["ctrl","expt"]; % want to plot both control and experiment in one plot
for q = 1:length(GLs) % shallow, deep
    % start extracting data
    for j = 1:n_simu
        % read the model
        group = folder_dir_groups{q};
        % load both the model data and extracted centerline data
        dH_two = cell(1,2); 
        for tp = 1:length(md_types) % contrl and experiment
            md_type = md_types(tp);
            switch md_type
                case 'ctrl'
                    load([group.folder{j},'/', group.name{j}, '/', ctrl_name])
                case 'expt' % here we use effective pressure experiment by default
                    load([group.folder{j},'/', group.name{j}, '/', expt_name])
                otherwise
                    warning('unsupported input')
            end
            % find the average thinning rate along thalweg
            idx_first = 1; idx_last = size(md.results.TransientSolution,2);
            [H_first,x,y] = plot_thalweg_profile(md, ds, idx_first, 'thickness','ice',idx_last);
            [H_last,~,~]  = plot_thalweg_profile(md, ds, idx_last,  'thickness','ice',idx_last);
            dH_two{tp} = H_first - H_last;
        end
        modelname = md.miscellaneous.name(9:end);
        [Ws_a(q,j), GLs_a(q,j), FCs_a(q,j)] = parse_modelname(modelname);
        
        dH = dH_two{2} - dH_two{1}; % isolate dH from effective pressure
        maxdHdt(q,j) = max(dH)/t_total;
        % find attenuation distance
        cdf_dist = cumsum(dH)/nansum(dH);
        cdf_dist(isnan(cdf_dist)) = 1.0;
        [cdf_dist_u, u_i] = unique(cdf_dist); % keep only unique values for interpolation later
        x_u = x(u_i); %...and their corresponding x values
        dL(q,j) = sim_params.terminus0_x - interp1(cdf_dist_u,x_u,1/exp(1)); % attenuation dist. wrt. to initial calving front position

        % print
        disp(['Model: ' modelname])
        disp(['Max dh/dt maximum is ' num2str(maxdHdt(q,j)) ' m/a'])
        disp(['Attenuation distance is ' num2str(dL(q,j)) ' m'])
    end
end

%% Make a grid plot
% Marker is filled circle, where 
% ... the size represents the max dh/dt
% ... the color represents the attenuation distance

grid_w = [5000,5000,5000;8000,8000,8000;11000,11000,11000];
grid_fc = [3e4, 6e4, 12e4;3e4, 6e4, 12e4;3e4,  6e4,  12e4];
[grid_X, grid_Y] = meshgrid(1:3,3:-1:1);
% marker style specs
size_range = [5,9];
red = [255, 51, 153]/255;
magenta = [153, 153, 255]/255;
color_range = [red; magenta];
maxdHdt_range = [min(maxdHdt,[],'all'), max(maxdHdt,[],'all')];
dL_range = [min(dL,[],'all'), max(dL,[],'all')];

% __ glacier with floating terminus: choose your depth
% two options: 0 meter or 400 meter
depth = 400;
figure('Position',[100,100,600,600])
for i = 1:length(grid_w(:))
    disp(['W = ' num2str(grid_w(i)) ', FC = ' num2str(grid_fc(i)) ' , GL = ' num2str(depth)])
    w_bl = grid_w(i) == Ws_a; 
    fc_bl = grid_fc(i) == FCs_a;
    gl_bl = depth == GLs_a;
    idx = find(w_bl + fc_bl + gl_bl == 3);

    % marker spec
    color_interp = interp1(dL_range, color_range, dL(idx));
    size_interp = interp1(maxdHdt_range, size_range, maxdHdt(idx));
    % coordinate
    scatter(grid_X(i), grid_Y(i), exp(size_interp), color_interp, 'filled');
    hold on
end
xlim([-1,5]); ylim([-1,5]);
plotname = ['scatter_effP_' num2str(depth) '.pdf'];

exportgraphics(gcf,['plots/composite_effP/' plotname],'ContentType','vector')