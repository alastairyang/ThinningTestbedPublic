% plot the scatter plot: lagtime gradient and lagtime magnitude
% first plot is calving-only
reversal_year = 5+8; % 5 year stationary in the beginning and 8 accelerated retreat years.
sample_number = 20; % make # samples along the center line
sample_interval = 1000; % meter; distance between two control points

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
        scp1 = scatter(slope, mean_tau, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        %errorbar(slope, mean_tau, min_tau, max_tau,'k')
        hold on

    end
    ylim([0,12])
    xlim([0,0.6])
    ylabel('$\tilde{\tau}$ (yr)','Interpreter','latex','FontSize',16)
    xlabel('$\partial \tau/\partial x$ (yr/km)','Interpreter','latex','FontSize',16)

end
legend off

% save to folder
saveas(gcf, 'plots/calve_timelag_scatter.pdf')

% set all the data in the previous plot more transparent
alpha(0.2)

% make calving-only less transparent; add mass unloading + calving

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

%figure('Position',[100,100,500,500]);
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
        %errorbar(slope, mean_tau, min_tau, max_tau,'k')
        hold on

    end
    ylim([0,12])
%     subplot(3,3,8)
%     xlabel('Distance to ice front (km)','Interpreter','latex',FontSize=13)
%     subplot(3,3,4)
%     ylabel('Lag time (yr)','Interpreter','latex',FontSize=13)
end
legend off

% save to folder
saveas(gcf, 'plots/mu_calve_timelag_scatter.pdf')


%% functions
function bool = compare_GLvalue(str, val)
    filename_split = split(str, '_');
    initials = string(cellfun(@(s) s(1:2), filename_split, 'UniformOutput', false));
    GLvalue = filename_split(strcmp('GL', initials));
    GLvalue = str2double(GLvalue{1}(3:end));

    bool = GLvalue == val;
end