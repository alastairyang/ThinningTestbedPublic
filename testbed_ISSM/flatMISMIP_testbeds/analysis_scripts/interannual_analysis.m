%% Analyzing timeseries from inter-annual basal perturbation experiment
% Here analyze the h(t) from inter-annual basal gaussian perburbation. I
% detrend each time series with polynomial fitting. For the cyclic
% component, I find the amplitude of the oscillation and characterize how
% this changes along the glacier flow line as an "attenuation length
% scale". On the other hand, the local perturbation can cause changes in
% the mean state (trend), which we also quantify.

% Here we use the data sampled and saved from the model class (saved in
% "analyzed_data" folder). We look at the experiment and control runs
% simultaneously.
perturb_index = 2;
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

ctrl_foldername = 'analyzed_data/calve_only';
expt_foldername = 'analyzed_data/gp1_calve';
ctrl_folder_prefix = 'ht_calve_';
expt_folder_prefix = 'ht_gp1_';
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

% shallower grounding line
n_simu = size(expt_folder_dir_groups{shallowGL_i}, 1);
% initialize
A0s_shallow = zeros(1, n_simu);
decay_lengths_shallow = zeros(1,n_simu);
for j = 1:n_simu
    % load in data from paths
    ctrl_path = ctrl_folder_dir_groups{shallowGL_i}(j,:);
    expt_path = expt_folder_dir_groups{shallowGL_i}(j,:);
    ctrl = load([ctrl_path.folder{1},'/', ctrl_path.name{1}]);
    expt = load([expt_path.folder{1},'/', expt_path.name{1}]);
    % start analysis
    [A0s_shallow(j), decay_lengths_shallow(j)] = gauss_perturb_analysis(expt.ht_data.h, ctrl.ht_data.h, ctrl.ht_data.x, gauss_xloc);

    % assign plot legends according to parameters
    [W, GL, FC] = parse_modelname(ctrl_path.name{1});
    W_symb = Ws_symb(W==Ws); % marker size
    GL_symb = GLs_symb(GL==GLs); % marker type
    FC_symb = FCs_symb(FC==FCs,:); % marker color (rgb)
    scatter(A0s_shallow(j), decay_lengths_shallow(j), W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
    hold on;
end

%% Deeper grounding line
% deeper grounding line
n_simu = size(expt_folder_dir_groups{deeperGL_i}, 1);
% initialize
A0s_deep = zeros(1, n_simu);
decay_lengths_deep = zeros(1,n_simu);
figure;
for j = 1:n_simu
    % load in data from paths
    ctrl_path = ctrl_folder_dir_groups{deeperGL_i}(j,:);
    expt_path = expt_folder_dir_groups{deeperGL_i}(j,:);
    ctrl = load([ctrl_path.folder{1},'/', ctrl_path.name{1}]);
    expt = load([expt_path.folder{1},'/', expt_path.name{1}]);
    % start analysis
    [A0s_deep(j), decay_lengths_deep(j)] = gauss_perturb_analysis(expt.ht_data.h, ctrl.ht_data.h, ctrl.ht_data.x, gauss_xloc);

    % assign plot legends according to parameters
    [W, GL, FC] = parse_modelname(ctrl_path.name{1});
    W_symb = Ws_symb(W==Ws); % marker size
    GL_symb = GLs_symb(GL==GLs); % marker type
    FC_symb = FCs_symb(FC==FCs,:); % marker color (rgb)
    scatter(A0s_deep(j), decay_lengths_deep(j), W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
    hold on;
end

%% APPENDIX: functions
function [A0, decay_length] = gauss_perturb_analysis(perturb_ht, nopertb_ht, pos, gauss_xloc)
%GAUSS_PERTURB_ANALYSIS analyze the time series from inter-annual basal
%perturbation. 
%
%   Input:
%       perturb_ht [double, array]: h(t) for the perturbation run
%       nopertb_ht [double, array]: h(t) for the control run
%       pos        [double, array]: distance of the control point in meter from x=0 
%       gauss_xloc [double]: distance of the center of gaussian
%       perturbation from x = 0
%
%   Output:
%       A0           [double]: amplitude of the cyclic component
%       decay_length [double]: length scale of the attentuation in
%       amplitude

    n_gcp = length(pos);
    % Decompose the timeseries
    % find the difference between control and experiment
    delta_ht = perturb_ht - nopertb_ht;
    % crop out the period where there is no inter-annual perturbation
    yr = 5; dt = 0.1;
    t = yr/dt;
    % decompose timeseries into short term (cyclic) and long term (trend)
    delta_ht = delta_ht(t+1:end-t,:);
    STs = detrend(delta_ht, 3);
    LTs = delta_ht - STs;
    
    % make a plot of the decomposed trend and inter-annual variability
    fig = figure;
    colororder(cool(size(LTs,2)))
    subplot(1,2,1)
    plot(1:size(LTs,1), LTs);
    subplot(1,2,2)
    plot(1:size(STs,1), STs);
%     
    % find peaks and troughs
    peaks = zeros(n_gcp,1);
    for i = 1:n_gcp
        pks = findpeaks(STs(:,i));
        troughs = findpeaks(-1*STs(:,i));
        peaks(i) = mean(pks) + mean(troughs);
    end
    
    % find the control point immmediately downstream of the center of
    % the gaussian perturbation patch
    perturb_xi = find(pos > gauss_xloc,1,'last') - 3;
    peaks_downstream = peaks(1:perturb_xi);
    % estimate the decay length scale
    %options = optimset('PlotFcns',@optimplotfval);
    options = optimset('Display','notify');
    decay_v = [0.00001]; % decay length scale initial guess
    dist_to_perturb = pos(1:perturb_xi) - pos(perturb_xi);
    A0 = max(peaks_downstream);
    decay_length = fminsearch(@minimize_attenuation, decay_v, options, A0, dist_to_perturb, peaks_downstream);
    % plot the attenuation model and the data
    %plot_attenuation(decay_length, A0, dist_to_perturb, peaks_downstream)
end

function err = minimize_attenuation(v, A, x, data)
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);
    err = sqrt(mean((A_decay' - data).^2));
end

function plot_attenuation(v, A, x, data)
    beta = v(1); % decay length scale
    A_decay = A.*exp(-beta.*x);

    figure;
    scatter(x, data, 4,'red'); hold on;
    plot(x, A_decay, '-b'); hold off
    legend(["data","model"])

end


