%% Explore the proper degree of polynomials
% When using the polynomials to decompose the time series, what is the best
% degree to use?
%% Mainscript
perturb_index = 8;
n_gcp = 20;
gauss_xloc = 3.2e4; % gaussian patch location from ice divide
glacier_length = 56500; % initial glacier length (x = 0 to calving front)
sample_interval = 1000; % meter; distance between two control points

% model parameters and plot parameters
% read in the model parameter table
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));

% specify the scatter plot data symbols
% GL: circle vs dot; FC: color, light to dark; W: size of the symbol
Ws_symb = [40,80,110];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

ctrl_foldername = 'analyzed_data/mu_calve';
expt_foldername = 'analyzed_data/mu_diffu_calve';
ctrl_folder_prefix = 'ht_mu_calve_';
expt_folder_prefix = 'ht_mu_diffu_';
ctrl_folder_dir = natsortfiles(dir([pwd '/' ctrl_foldername]));
expt_folder_dir = natsortfiles(dir([pwd '/' expt_foldername]));
ctrl_folder_dir = struct2table(ctrl_folder_dir);
expt_folder_dir = struct2table(expt_folder_dir);
% remove  '.' and '..'
bools = cellfun(@(s) ~strcmp(s(1),'.'), ctrl_folder_dir.name);
ctrl_folder_dir = ctrl_folder_dir(bools,:);
bools = cellfun(@(s) ~strcmp(s(1),'.'), expt_folder_dir.name);
expt_folder_dir = expt_folder_dir(bools,:);

plot_idx = 0;

ylabel_i = [1,4,7];
xlabel_i = [7,8,9];

% split the folder_dir into two groups, separated by grounding line depth
% control run
ctrl_folder_dir_groups = cell(1,2);
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(ctrl_folder_dir,1),1);
    for j = 1:size(ctrl_folder_dir.name)
        GL_bool(j) = compare_GLvalue(ctrl_folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    ctrl_folder_dir_groups{i} = ctrl_folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end
% experiment run
expt_folder_dir_groups = cell(1,2);
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(expt_folder_dir,1),1);
    for j = 1:size(expt_folder_dir.name)
        GL_bool(j) = compare_GLvalue(expt_folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    expt_folder_dir_groups{i} = expt_folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end


%% Time series analysis
% we divide the dicussions by the grounding line depth
[~, shallowGL_i] = min(GLs);
[~, deeperGL_i]  = max(GLs);

%% Analyze
geom_type = 'deep';
degrees = 4:2:18;

switch geom_type
    case 'deep'
        geom_i = deeperGL_i;
    case 'shallow'
        geom_i = shallowGL_i;
    otherwise
        warning('unknown geom type')
end

n_simu = size(expt_folder_dir_groups{geom_i}, 1);
% initialize
A0s_shallow = zeros(1, n_simu);
decay_lengths_shallow = zeros(1,n_simu);
freqs = cell(1, n_simu);
rel_t = cell(1, n_simu);
x_pos = cell(1, n_simu);
STs_errors = zeros(n_gcp, length(degrees));
LTs_errors = zeros(n_gcp, length(degrees));
simu_errors = cell(2,n_simu); % save both errors for all simulations
for j = 1:n_simu % iterate over simulation for the given depth
    for i_n = 1:length(degrees) % iterate over polynomial degree
        % load in data from paths
        n = degrees(i_n);
        ctrl_path = ctrl_folder_dir_groups{geom_i}(j,:);
        expt_path = expt_folder_dir_groups{geom_i}(j,:);
        ctrl = load([ctrl_path.folder{1},'/', ctrl_path.name{1}]);
        expt = load([expt_path.folder{1},'/', expt_path.name{1}]);
        % get the last grounding line position
        % we only look at the data points beteen that and the localized
        % gaussian perturbation
        gl_last_xloc = min(ctrl.ht_data.gl(end), expt.ht_data.gl(end));
        % start analysis
        [STs, LTs] = gauss_perturb_analysis(expt.ht_data.h, ctrl.ht_data.h, n);
        STs_errors(:,i_n) = transpose(mean(abs(STs(1:50,:) - 0),1)); % the mean is taken for the first 50 timesteps
        LTs_errors(:,i_n) = transpose(mean(abs(LTs(1:50,:) - 0),1)); 
        disp(['model no.',num2str(j), ' is completed'])
    end
    simu_errors{1,j} = STs_errors;
    simu_errors{2,j} = LTs_errors;
end

%%
figure
for k = 1:n_simu
    subplot(2,1,1); 
    fanChart(degrees, transpose(simu_errors{1,k}), 'mean',[10,50,90], 'alpha',0.2,'colormap',{'shadesOfColor',[0 0 .8]}); hold on
    subplot(2,1,2); 
    fanChart(degrees, transpose(simu_errors{2,k}), 'mean',[10,50,90], 'alpha',0.2,'colormap',{'shadesOfColor',[0 0 .8]}); hold on
end
subplot(2,1,1); title('Error in cyclic component','Interpreter','latex','FontSize',13)
subplot(2,1,2); title('Error in trend component','Interpreter','latex','FontSize',13);
xlabel('degree n polynomial','Interpreter','latex','FontSize',13)
%exportgraphics(gcf,'plots/poly_error.pdf','ContentType','vector')

%% APPENDIX: Functions
function [STs, LTs] = gauss_perturb_analysis(perturb_ht, nopertb_ht, n)
%GAUSS_PERTURB_ANALYSIS analyze the time series from inter-annual basal
%perturbation. 
%
%   Input:
%       perturb_ht [double, array]: h(t) for the perturbation run
%       nopertb_ht [double, array]: h(t) for the control run
%       n [int]: the degree of polynomials


    % Decompose the timeseries
    % find the difference between control and experiment
    if any(size(nopertb_ht) ~= size(perturb_ht))
        if size(nopertb_ht, 1) == size(perturb_ht,1) + 1
            % we allow for one step difference; if true, crop one
            nopertb_ht = nopertb_ht(1:end-1,:);
        else
            error('Sampled dh/dt time dimension incompatible!')
        end
    end
    delta_ht = perturb_ht - nopertb_ht;
    STs = detrend(delta_ht, n);
    LTs = delta_ht - STs;
end