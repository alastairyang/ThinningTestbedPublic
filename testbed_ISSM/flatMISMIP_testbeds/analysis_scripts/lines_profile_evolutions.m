%% Plot the glacier lateral profiles over time
md1 = load('long_models_yang/model_W11000_GL400_FC30000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;
md2 = load('long_models_yang/model_W5000_GL400_FC120000/MISMIP_yangTransient_Calving_MassUnloading.mat').md;

%% Recover the glacier profiles overtime
dt = 0.1;
t_interval = 2;
ds = 50;
[surface1, base1, bed1, profiles_t1] = sample_profile_evol(md1, dt, ds, t_interval);
[surface2, base2, bed2, profiles_t2] = sample_profile_evol(md2, dt, ds, t_interval);

%% plot
% load colormap
cp = load('plots/colormap/lajolla.mat').lajolla;

% make the plot
plot_x = 0:ds:(size(surface1,2)-1)*ds;
figure('Position',[100,100,900,600])
tiledlayout(2,1,"TileSpacing","none")
nexttile
plot(plot_x/1000, surface1,'LineWidth',1); hold on;
plot(plot_x/1000, base1,'LineWidth',1); hold on;
plot(plot_x/1000, bed1,'-k','LineWidth',2); hold off
ylim([-600,1500])
line_colors = colormap_to_colororder(cp, size(surface1,1),1,100);
colororder(line_colors)
set(gca,'Xtick',[])
ylabel('Elevation (m)','FontName','Aria','FontSize',15)
text(45,1380,'Width = 11km, k = 3e4','FontName','Aria','FontSize',13)

nexttile
plot(plot_x/1000, surface2,'LineWidth',1); hold on;
plot(plot_x/1000, base2,'LineWidth',1); hold on;
plot(plot_x/1000, bed2,'-k','LineWidth',2); hold off
ylim([-600,1500])
line_colors = colormap_to_colororder(cp, size(surface2,1),1,100);
colororder(line_colors)

xlabel('Along-flow distance (km)','FontName','Aria','FontSize',15)
ylabel('Elevation (m)','FontName','Aria','FontSize',15)
text(45,1380,'Width = 5km, k = 12e4','FontName','Aria','FontSize',13)

exportgraphics(gcf,'plots/profiles_evolve.png','Resolution',600)

