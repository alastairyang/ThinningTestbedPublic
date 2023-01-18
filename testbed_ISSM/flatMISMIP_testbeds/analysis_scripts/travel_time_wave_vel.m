%% Wave travel time along the center flowline and average propagation velocity
% In this script we use the sampled dh/dt along the center flow line to
% derive
%   1. the travel time of the wave based off the forcing extrema
%   2. the average wave velocity; performing a linear regression on 1)

%%
% first plot is calving-only
reversal_year = 5+8; % 5 year stationary in the beginning and 8 accelerated retreat years.
sample_number = 20; % make # samples along the center line
sample_interval = 1000; % meter; distance between two control points

% read in model parameters
md_vars = readtable('md_var_combinations.csv');
Ws = sort(unique(md_vars.('fjord_width')));
GLs = sort(unique(md_vars.('delta_groundingline_depth')));
FCs = sort(unique(md_vars.('background_friccoef')));

% data directory
ctrl_foldername = 'analyzed_data/calve_only';
expt_foldername = 'analyzed_data/mu_calve';
ctrl_folder_prefix = 'ht_calve_';
expt_folder_prefix = 'ht_mu_calve_';
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

% specify the scatter plot data symbols
% GL: circle vs dot; FC: color, light to dark; W: size of the symbol
Ws_symb = [40,80,110];
GLs_symb = ["square","o"];
FCs_symb = [166,32,232;232,32,199;232,32,72]/255;

% preallocate
expt_folder_dir_groups = cell(1,2);
ctrl_folder_dir_groups = cell(1,2);

% control
% split the folder_dir into two groups, separated by grounding line depth
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(ctrl_folder_dir,1),1);
    for j = 1:size(ctrl_folder_dir.name)
        GL_bool(j) = compare_GLvalue(ctrl_folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    ctrl_folder_dir_groups{i} = ctrl_folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end
% experiment
% split the folder_dir into two groups, separated by grounding line depth
for i = 1:length(GLs)
    % skip the irrelevant ones
    GL_bool = zeros(size(expt_folder_dir,1),1);
    for j = 1:size(expt_folder_dir.name)
        GL_bool(j) = compare_GLvalue(expt_folder_dir.name(j), GLs(i));
    end
    % save the respective folder items to a cell
    expt_folder_dir_groups{i} = expt_folder_dir(find(GL_bool),:); %#ok<FNDSB> 

end

%% For control group
h = tiledlayout(2,1, 'TileSpacing','none','Padding','none');
%figure('Position',[100,100,500,500]);
for j = 1:length(GLs) % iterate over grounding line depths
    ctrl_folder_dir = ctrl_folder_dir_groups{j};
    md_count = 0;

    nexttile
    for i = 1:size(ctrl_folder_dir,1)
        md_count = md_count + 1;

        md = load(string(ctrl_folder_dir.folder(i))+"/"+ string(ctrl_folder_dir.name(i))).ht_data;
        modelname = ctrl_folder_dir.name(i);
        modelname = modelname{1}(length(ctrl_folder_prefix)+1:end-4);
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
        travel_time = min_dhdt_year;
        % create along flowline sampling distance
        % here we define the distance as to the first kept sampled point
        % (closest to the terminus)
        distances = (md.x(1) - md.x)/1000; 

        % skip the control points that are on the floating section at the
        % end of the simulation
        distances = distances(~(md.x > (md.gl(end)-1500)));
        travel_time = travel_time(~(md.x > (md.gl(end)-1500)));
        travel_time = travel_time - travel_time(1);
%         if j == length(GLs)
%             legend(["GL = "+string(num2str(GLs(1))), "GL = "+string(num2str(GLs(2)))]);
%         end
%         hold on
%         subplot_title = strrep(modelname(1:end-4), '_',', ');

        % more compact illustration: 1 scatter plot
        % get the linear fit to lag time - distance 
        P = polyfit(distances, travel_time, 1);
        slope = P(1); 
        intercpt = P(2);
        distances_forplot = distances(1):0.1:distances(end);
        travel_time_md = polyval(P, distances_forplot);
        % get the corresponding symbols for this scatter plot
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        % plot both the scatter and the linear fit
        scatter(distances, travel_time, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        hold on
        plot(distances_forplot, travel_time_md, '-.k'); 
        xlim([0,20])
        ylim([0,8])
        alpha(0.7)
    end
    ylabel(' Arrival time (yr)','Interpreter','latex','FontSize',16)
    xlabel('Distance to calving front $x$ (km)','Interpreter','latex','FontSize',16)
    % create fititious legends
    s1 = scatter(nan,nan, Ws_symb(1), 'k',GLs_symb(j)); hold on;
    s2 = scatter(nan,nan, Ws_symb(2), 'k',GLs_symb(j)); hold on;
    s3 = scatter(nan,nan, Ws_symb(3), 'k',GLs_symb(j)); hold on;
%     symb_names = {'narrow','standard','wide'};
%     [~, objh] = legend([s1 s2 s3], symb_names,'FontSize',12,'Box','off','Interpreter','latex','Location','northwest');
%     objhl = findobj(objh, 'type','scatter');
    %set(objhl, 'Markersize', Ws_symb(jj))
end
exportgraphics(gcf,'plots/ctrl_travel_time.pdf','ContentType','vector')
%% For experiment
h = tiledlayout(2,1, 'TileSpacing','none','Padding','none');
%figure('Position',[100,100,500,500]);
for j = 1:length(GLs) % iterate over grounding line depths
    expt_folder_dir = expt_folder_dir_groups{j};
    md_count = 0;

    nexttile
    for i = 1:size(expt_folder_dir,1)
        md_count = md_count + 1;

        md = load(string(expt_folder_dir.folder(i))+"/"+ string(expt_folder_dir.name(i))).ht_data;
        modelname = expt_folder_dir.name(i);
        modelname = modelname{1}(length(expt_folder_prefix)+1:end-4);
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
        travel_time = min_dhdt_year;
        % create along flowline sampling distance
        % here we define the distance as to the first kept sampled point
        % (closest to the terminus)
        distances = (md.x(1) - md.x)/1000; 

        % skip the control points that are on the floating section at the
        % end of the simulation
        distances = distances(~(md.x > (md.gl(end)-1500)));
        travel_time = travel_time(~(md.x > (md.gl(end)-1500)));
        travel_time = travel_time - travel_time(1);
%         if j == length(GLs)
%             legend(["GL = "+string(num2str(GLs(1))), "GL = "+string(num2str(GLs(2)))]);
%         end
%         hold on
%         subplot_title = strrep(modelname(1:end-4), '_',', ');

        % more compact illustration: 1 scatter plot
        % get the linear fit to lag time - distance 
        P = polyfit(distances, travel_time, 1);
        slope = P(1); 
        intercpt = P(2);
        distances_forplot = distances(1):0.1:distances(end);
        travel_time_md = polyval(P, distances_forplot);
        % get the corresponding symbols for this scatter plot
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        % plot both the scatter and the linear fit
        scatter(distances, travel_time, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        hold on
        plot(distances_forplot, travel_time_md, '-.k'); 
        xlim([0,20])
        ylim([0,8])
        alpha(0.7)
    end
    ylabel(' Arrival time (yr)','Interpreter','latex','FontSize',16)
    xlabel('Distance to calving front $x$ (km)','Interpreter','latex','FontSize',16)
    % create fititious legends
    s1 = scatter(nan,nan, Ws_symb(1), 'k',GLs_symb(j)); hold on;
    s2 = scatter(nan,nan, Ws_symb(2), 'k',GLs_symb(j)); hold on;
    s3 = scatter(nan,nan, Ws_symb(3), 'k',GLs_symb(j)); hold on;
%     symb_names = {'narrow','standard','wide'};
%     [~, objh] = legend([s1 s2 s3], symb_names,'FontSize',12,'Box','off','Interpreter','latex','Location','northwest');
%     objhl = findobj(objh, 'type','scatter');
    %set(objhl, 'Markersize', Ws_symb(jj))
end
exportgraphics(gcf,'plots/mu_travel_time.pdf','ContentType','vector')

%% save to folder
saveas(gcf, 'plots/calve_timelag_scatter.pdf')

% set all the data in the previous plot more transparent
