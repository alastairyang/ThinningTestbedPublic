%% dhdt_gifs.m
% in this script, I plot GIFs of retreating glaciers (its lateral
% profile) along with the dh/dt at the selected points. I compare the
% glacier with retreat-only and with retreat + effective pressure feedback

%% Parameter
sample_interval = 3000; % meter; distance between two control points
sample_number = 5; % make # samples along the center line
ds = 50; % grid spacing for mesh->grid interpolation, meter
foldername = 'long_models_yang';

calving_modelname = 'MISMIP_yangTransient_CalvingOnly.mat';
calving_mu_modelname = 'MISMIP_yangTransient_Calving_MassUnloading.mat';

modelname = 'model_W11000_GL400_FC120000';

%% load data
md1_path = string([foldername,'/', modelname,'/',calving_modelname]);
md2_path = string([foldername,'/', modelname,'/',calving_mu_modelname]);
md1 = load(md1_path).md;
md2 = load(md2_path).md;
mds = [md1, md2];

% Extract data
% we need to extract profile(t) and dh/dt(t) from ISSM model class. While
% there are saved dh/dt data in "analyzed_data" folder, to make this script
% standalone we are not using them.

ds_i = sample_interval/ds;
% store in a cell array
mds_data = cell(1,2);
for i = 1:length(mds)
    md = mds(i);
    nt = size(md.results.TransientSolution,2);

    % geometry parameters
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
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
    mesh_elements = md.mesh.elements;
    mesh_x = md.mesh.x;
    mesh_y = md.mesh.y;
    % retrieve the last position of the ice front; spclevelset ~= 0
    pos = find(md.levelset.spclevelset(1:end-1,end) < 1 & md.levelset.spclevelset(1:end-1,end) > -1);
    x_front = min(mesh_x(pos));

    % find grid index where we want to sample h(t)
    [~, x_i_nearest] = min(abs(x - x_front));
    % sampled points: start at __ behind the last ice front
    front_i = x_i_nearest-ds_i;
    end_i   = front_i - sample_number*ds_i;
    sample_i = front_i:-ds_i:(end_i+ds_i);

    % remove model class; data store in table instead to clear space
    results = struct2table(md.results.TransientSolution);
    empty_md = md;
    empty_md.results.TransientSolution = [];

    % initialize space to store h(t)
    thalweg_sample_ht = [];
    surface_t = [];
    base_t = [];

    % get bed data
    bed = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);
    % get all dh/dt(t), surface(t), base(t) data
    for j = 1:nt
        surface = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
            results.Surface{j},...
            x, y, NaN);
        base    = InterpFromMeshToGrid(empty_md.mesh.elements, mesh_x, mesh_y,...
            results.Base{j},...
            x, y, NaN);
        surface_points  = surface(mid_i,sample_i);
        surface_profile = surface(mid_i,:);
        base_profile     = base(mid_i,:);
        thalweg_sample_ht = [thalweg_sample_ht; surface_points];
        surface_t = [surface_t; surface_profile];
        base_t = [base_t; base_profile];
    end
    % time vector
    time = results.time;
    % shift time to start at 0
    time = time(1:end-1);
    time = time - time(1);
    dt = mean(time(2:end) - time(1:end-1));

    thalweg_sample_ht = thalweg_sample_ht - thalweg_sample_ht(1,:);
    % crop the one extra time step in simulation
    thalweg_sample_ht = thalweg_sample_ht(1:end-1,:);
    thalweg_sample_dhdt = (thalweg_sample_ht(2:end,:) - thalweg_sample_ht(1:end-1,:))/dt;
    % smooth the dh/dt data
    for k = 1:size(thalweg_sample_dhdt,2)
        thalweg_sample_dhdt(:,k) = smooth(thalweg_sample_dhdt(:,k),20);
    end

    % organize into data structure
    data.h = thalweg_sample_dhdt;
    data.t = time;
    data.bed = bed_profile;
    data.base = base_t;
    data.surface = surface_t;

    % save to the cell array
    mds_data{i} = data;
end
clear mds
%% Plot calving only
% data
data = mds_data{1};
gif_name = 'calving_gif.gif';
dhdt_ylim = -12;

% make the along flow x-axis for profile plots
x = 0:ds:ds*length(bed_profile)-1;
% dhdt coloring
dhdt_color_length = sample_number;
lightblue = [80, 203, 253]/255;
darkblue = [0, 128, 255]/255;
dhdt_colors_p = [linspace(lightblue(1),darkblue(1),dhdt_color_length)',...
                linspace(lightblue(2),darkblue(2),dhdt_color_length)',...
                linspace(lightblue(3),darkblue(3),dhdt_color_length)'];

% profile coloring
profile_color_length = 260;
grey = [96, 96, 96]/255;
darkgrey = [32, 32, 32]/255;
profile_colors_p = [linspace(grey(1),darkgrey(1),profile_color_length)',...
                    linspace(grey(2),darkgrey(2),profile_color_length)',...
                    linspace(grey(3),darkgrey(3),profile_color_length)'];
% semi-fake retreat sequence
retreat_sequence = [0,0,0,0,0,0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0,0,0,0,0]*200;
fake_t = 0:length(retreat_sequence)-1;
fake_t_interp = 0:0.1:length(retreat_sequence)-0.1;
retreat_seq_interp = interp1(fake_t, retreat_sequence, fake_t_interp);

% parameters
start_i = 50;
end_i = length(data.t)-1;
gcp_x = ds*sample_i;
profile_y_upperlimit = 2000;

figure;
% we first adjust the positions of the two subplots
% then call their handle to plot. The reason: using "set" command can clear
% the the "hold on" command.
p1 = subplot(2,1,1);
p2 = subplot(2,1,2);
q1 = get(p1, 'Position');
q2 = get(p2, 'Position');
q1(2) = q2(2)+q2(4);
set(p1, 'pos', q1);

gif(['plots/',gif_name],'DelayTime',0.4)
for i = start_i:10:end_i
    %subplot(2,1,1);
    % plot the bed
    axes(p1); % first subplot
    plot(x, data.bed,'k'); hold on;
    % plot the evolving profile
    plot(x, data.surface(i,:),'Color',profile_colors_p(i,:)); hold on
    plot(x, data.base(i,:),'Color',profile_colors_p(i,:)); hold on;
    % plot the x locations of the control points
    for k = 1:length(gcp_x)
        scatter(gcp_x(k), profile_y_upperlimit,55,dhdt_colors_p(k,:),'filled'); hold on;
    end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
    % plot the retreat sequence
    axis off
    axes('Position',[0.8,0.6,0.1,0.1])
    ylim([0,10])
    plot(fake_t_interp(start_i:i), retreat_seq_interp(start_i:i),'Color',[255,128,0]/255,'LineWidth',2)
    xlim([0,26])
    box on
    axis off
    ylim([-600,profile_y_upperlimit])

    % plot the dh/dt(t)
    %subplot(2,1,2);
    axes(p2)
    h = plot(data.t(start_i:i), transpose(data.h(start_i:i,:)),'LineWidth',1.5); hold on
    set(h, {'color'}, num2cell(dhdt_colors_p, 2));
    ylim([dhdt_ylim,0])
    xlim([data.t(start_i),data.t(end_i)])
    xlabel('Time (yr)','Interpreter','latex','FontSize',16)
    ylabel('dH/dt (m/a)','Interpreter','latex','FontSize',16)

    gif
    pause(0.1)

end

%% Plot calving + effective pressure feedback
% data
data = mds_data{2};
gif_name = 'mu_calving_gif.gif';
dhdt_ylim = -55;

% make the along flow x-axis for profile plots
x = 0:ds:ds*length(bed_profile)-1;
% dhdt coloring
dhdt_color_length = sample_number;
lightblue = [80, 203, 253]/255;
darkblue = [0, 128, 255]/255;
dhdt_colors_p = [linspace(lightblue(1),darkblue(1),dhdt_color_length)',...
                linspace(lightblue(2),darkblue(2),dhdt_color_length)',...
                linspace(lightblue(3),darkblue(3),dhdt_color_length)'];

% profile coloring
profile_color_length = 260;
grey = [96, 96, 96]/255;
darkgrey = [32, 32, 32]/255;
profile_colors_p = [linspace(grey(1),darkgrey(1),profile_color_length)',...
                    linspace(grey(2),darkgrey(2),profile_color_length)',...
                    linspace(grey(3),darkgrey(3),profile_color_length)'];
% semi-fake retreat sequence
retreat_sequence = [0,0,0,0,0,0,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,0,0,0,0,0]*200;
fake_t = 0:length(retreat_sequence)-1;
fake_t_interp = 0:0.1:length(retreat_sequence)-0.1;
retreat_seq_interp = interp1(fake_t, retreat_sequence, fake_t_interp);

% parameters
start_i = 50;
end_i = length(data.t)-1;
gcp_x = ds*sample_i;
profile_y_upperlimit = 2000;

figure;
% we first adjust the positions of the two subplots
% then call their handle to plot. The reason: using "set" command can clear
% the the "hold on" command.
p1 = subplot(2,1,1);
p2 = subplot(2,1,2);
q1 = get(p1, 'Position');
q2 = get(p2, 'Position');
q1(2) = q2(2)+q2(4);
set(p1, 'pos', q1);

gif(['plots/',gif_name],'DelayTime',0.4)
for i = start_i:10:end_i
    %subplot(2,1,1);
    % plot the bed
    axes(p1); % first subplot
    plot(x, data.bed,'k'); hold on;
    % plot the evolving profile
    plot(x, data.surface(i,:),'Color',profile_colors_p(i,:)); hold on
    plot(x, data.base(i,:),'Color',profile_colors_p(i,:)); hold on;
    % plot the x locations of the control points
    for k = 1:length(gcp_x)
        scatter(gcp_x(k), profile_y_upperlimit,55,dhdt_colors_p(k,:),'filled'); hold on;
    end
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
    % plot the retreat sequence
    axis off
    axes('Position',[0.8,0.6,0.1,0.1])
    ylim([0,10])
    plot(fake_t_interp(start_i:i), retreat_seq_interp(start_i:i),'Color',[255,128,0]/255,'LineWidth',2)
    xlim([0,26])
    box on
    axis off
    ylim([-600,profile_y_upperlimit])

    % plot the dh/dt(t)
    %subplot(2,1,2);
    axes(p2)
    h = plot(data.t(start_i:i), transpose(data.h(start_i:i,:)),'LineWidth',1.5); hold on
    set(h, {'color'}, num2cell(dhdt_colors_p, 2));
    ylim([dhdt_ylim,0])
    xlim([data.t(start_i),data.t(end_i)])
    xlabel('Time (yr)','Interpreter','latex','FontSize',16)
    ylabel('dH/dt (m/a)','Interpreter','latex','FontSize',16)

    gif
    pause(0.1)

end

