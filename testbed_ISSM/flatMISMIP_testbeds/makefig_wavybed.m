% This creates a .png file showing the rough bed geometry, colored by
% deviation of height from mean height.
load('random_beds/dD80_H07_RO4000_3.mat')
load('plots/colormap/vik.mat')
%% plot
figure('Position',[100,100,600,400])
f = mesh(rand_bed.X, rand_bed.Y, rand_bed.z);
view([320,50])
%f.FaceColor = 'flat';
clim([-250,250])
colormap(vik); cb = colorbar;
cb.Ticks = [-250,-150,-50,50,150,250];
cb.Label.String = "Elevation offset (m)";
cb.FontSize = 16; cb.Label.FontSize = 18;
axis off;

% export
exportgraphics(gcf,'plots/random_bed_mesh.png','Resolution',600)