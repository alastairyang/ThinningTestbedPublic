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
% pre-allocate
ctrl_slopes = zeros(length(GLs),size(ctrl_folder_dir,1));
ctrl_slopes_ci = zeros(length(GLs),size(ctrl_folder_dir,1),2);

figure('Position',[100,100,600,700]);
tiledlayout(2,1, 'TileSpacing','none','Padding','none');
titles = ["Shallow","Deep"];
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
        smooth_window = 30;
        % to smooth data, we have to iterate over each flow line
        % reshape the matrix into a vector and smooth will create
        % jumps/shocks
        dhdt = zeros(size(md.h,1)-1, size(md.h,2));
        for line_i = 1:size(md.h,2)
            md.h(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end

        dhdt_ori = diff(md.h,1,1)/(md.t(2) - md.t(1));
        for line_i = 1:size(dhdt,2)
            dhdt(:,line_i) = smooth(dhdt_ori(:,line_i), smooth_window);
        end
        

        % first interp and then find the lag year to calving front reversal
        t_finer = md.t(1):0.01:md.t(end);
        dhdt = interp1(md.t(1:end-1),dhdt,t_finer,'spline');
        [min_dhdt, min_dhdt_i] = min(dhdt,[], 1);
        min_dhdt_year = t_finer(min_dhdt_i);
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
        %[pp, ~] = regress(travel_time(:), [distances(:) ones(numel(distances),1)]);
        %distances_forplot = distances(1):0.1:distances(end);
        %travel_time_md = polyval(pp, distances_forplot);
        [pp, ~] = regress(distances(:), [travel_time(:) ones(numel(distances),1)]);
        travel_time_forplot = travel_time(1):0.01:travel_time(end);
        distances_forplot = polyval(pp, travel_time_forplot);
        % get the corresponding symbols for this scatter plot
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        % plot both the scatter and the linear fit
        %scatter(distances, travel_time, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        scatter(travel_time, distances, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        hold on
        %plot(distances_forplot, travel_time_md, '-.k'); 
        plot(travel_time_forplot, distances_forplot, '-.k')
        xlim([0,6])
        ylim([0,25])
        alpha(0.7)
        % linear regression again, but x y reversed
        [pp, pp_ci] = regress(distances(:), [travel_time(:) ones(numel(distances),1)]);
        slope = pp(1); 
        ctrl_slopes(j,i) = slope;
        ctrl_slopes_ci(j,i,:) = pp_ci(1,:);
    end
    xlabel(' Arrival time (yr)','Interpreter','latex','FontSize',16)
    ylabel('Distance to calving front (km)','Interpreter','latex','FontSize',16)
    text(5.5, 18, titles(j),'Interpreter','latex','FontSize',15)
end
exportgraphics(gcf,'plots/ctrl_travel_time.pdf','ContentType','vector')
%% For experiment
expt_slopes = zeros(length(GLs),size(expt_folder_dir,1));
expt_slopes_ci = zeros(length(GLs),size(expt_folder_dir,1),2);

figure('Position',[100,100,600,700])
tiledlayout(2,1, 'TileSpacing','none','Padding','none');
titles = ["Shallow","Deep"];

hCopys = []; % storing the scatter elements
kin_wave_vels = zeros(length(GLs)*size(expt_folder_dir,1), 4); % 4 columns: width, depth, coef, slope
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
        smooth_window = 30;
        % to smooth data, we have to iterate over each flow line
        % reshape the matrix into a vector and smooth will create
        % jumps/shocks
        dhdt = zeros(size(md.h,1)-1, size(md.h,2));
        for line_i = 1:size(md.h,2)
            md.h(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end

        dhdt_ori = diff(md.h,1,1)/(md.t(2) - md.t(1));
        for line_i = 1:size(dhdt,2)
            dhdt(:,line_i) = smooth(dhdt_ori(:,line_i), smooth_window);
        end

        % first interp and then find the lag year to calving front reversal
        t_finer = md.t(1):0.01:md.t(end);
        dhdt = interp1(md.t(1:end-1),dhdt,t_finer,'spline');
        [min_dhdt, min_dhdt_i] = min(dhdt,[], 1);
        min_dhdt_year = t_finer(min_dhdt_i);
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
        [pp, ~] = regress(distances(:), [travel_time(:) ones(numel(distances),1)]);
        travel_time_forplot = travel_time(1):0.01:travel_time(end);
        distances_forplot = polyval(pp, travel_time_forplot);
        % get the corresponding symbols for this scatter plot
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        % plot both the scatter and the linear fit
        %scatter(distances, travel_time, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        scatter(travel_time, distances, W_symb, GL_symb, 'MarkerFaceColor',FC_symb,'MarkerEdgeColor','k');
        hold on
        %plot(distances_forplot, travel_time_md, '-.k'); 
        plot(travel_time_forplot, distances_forplot, '-.k')
        xlim([0,6])
        ylim([0,25])
        alpha(0.7)
%         % add the legend with customized marker size 
%         hCopy = copyobj(h, ax); 
%         set(hCopy,'XData', NaN', 'YData', NaN);
%         hCopy.MarkerSize = W_symb*0.15;
%         hCopys = [hCopys, hCopy];

        % save to workspace
        % linear regression again, but x y reversed
        [pp, pp_ci] = regress(distances(:), [travel_time(:) ones(numel(distances),1)]);
        slope = pp(1); 
        expt_slopes(j,i) = slope;
        expt_slopes_ci(j,i,:) = pp_ci(1,:);

        % save the slopes (wave velocity) along with model parameters
        % to a table
        % first convert subscript indexing to linear indexing
        id = sub2ind([length(GLs), size(expt_folder_dir,1)], j, i);
        kin_wave_vels(id,1) = W;
        kin_wave_vels(id,2) = GL;
        kin_wave_vels(id,3) = FC;
        kin_wave_vels(id,4) = slope;

    end
    xlabel(' Arrival time (yr)','Interpreter','latex','FontSize',16)
    ylabel('Distance to calving front $x$ (km)','Interpreter','latex','FontSize',16)
    text(5.5, 18, titles(j),'Interpreter','latex','FontSize',15)

end
exportgraphics(gcf,'plots/mu_travel_time.pdf','ContentType','vector')

% make the kinematic wave estimate a table and export to a spreedsheet
kin_wave_vels = array2table(kin_wave_vels, 'VariableNames',["Width","Depth","Coefficient","KW Velocity (km/a)"]);
writetable(kin_wave_vels,'result_tables/kinematic_wave_estimate_control.csv')

%% Plot all the slopes
ctrl_x = 3;
expt_x = 5;
plot_xlim = [ctrl_x-1.5, expt_x+1.5];
xlabels = ["Shallow","Deep"];
figure('Position',[100,100,1000,700])
tiledlayout(1,2,"TileSpacing","none")
for j = 1:length(GLs)
    nexttile
    for i = 1:size(expt_folder_dir,1)
        % load the data to get model parameter information
        md = load(string(expt_folder_dir.folder(i))+"/"+ string(expt_folder_dir.name(i))).ht_data;
        modelname = expt_folder_dir.name(i);
        modelname = modelname{1}(length(expt_folder_prefix)+1:end-4);
        % get the model parameter
        [W, GL, FC] = parse_modelname(modelname);
        W_symb = Ws_symb(W==Ws);
        GL_symb = GLs_symb(GL==GLs);
        FC_symb = FCs_symb(FC==FCs,:);
        % add the data points
        % ctrl
        errorbar(ctrl_x, ctrl_slopes(j,i), squeeze(ctrl_slopes_ci(j,i,:))-ctrl_slopes(j,i), 'vertical', 'Color',FC_symb, 'Marker','.','MarkerSize',W_symb*0.4,'LineStyle','none');
        hold on;
        % expt
        errorbar(expt_x, expt_slopes(j,i), squeeze(expt_slopes_ci(j,i,:))-expt_slopes(j,i), 'vertical', 'Color',FC_symb, 'Marker','.','MarkerSize',W_symb*0.4,'LineStyle','none');
        hold on;
        % draw lines that connect
        plot([ctrl_x,expt_x],[ctrl_slopes(j,i),expt_slopes(j,i)],'Color',FC_symb,'LineStyle',':','LineWidth',W_symb*0.02);hold on
    end
    xlim(plot_xlim)
    ylim([0,20])
    xticks([ctrl_x, expt_x]);
    xticklabels({'without N','with N'});
    if j == 1; ylabel('Velocity (km/yr)','Interpreter','latex','FontSize',13); end
    if j == 2; set(gca, 'ytick',[]);end
    xlabel(xlabels(j),'Interpreter','latex','FontSize',13)
end
exportgraphics(gcf,'plots/wave_vels.png','Resolution',300)

%% save to folder
saveas(gcf, 'plots/calve_timelag_scatter.pdf')

%% 

% Create the plot
ax = axes(); 
hold on
h(1) = plot(linspace(1,5,25), rand(1,25), 'ro', 'DisplayName', 'foo');
h(2) = plot(1:5, rand(1,5), 'bo', 'DisplayName', 'bar');
% copy the objects
hCopy = copyobj(h, ax); 
% replace coordinates with NaN 
% Either all XData or all YData or both should be NaN.
set(hCopy(1),'XData', NaN', 'YData', NaN)
set(hCopy(2),'XData', NaN', 'YData', NaN)
% Note, these lines can be combined: set(hCopy,'XData', NaN', 'YData', NaN)
% To avoid "Data lengths must match" warning, assuming hCopy is a handle array, 
% use arrayfun(@(h)set(h,'XData',nan(size(h.XData))),hCopy)
% Alter the graphics properties
hCopy(1).MarkerSize = 30; 
hCopy(1).LineWidth = 2;
hCopy(2).MarkerSize = 3; 
% Create legend using copied objects
legend(hCopy)