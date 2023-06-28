% This creates .mat data file containing the randomly generated bed
% with fractal dimension of roughness specified by us.
deltaD = 80;
H = 0.7;
Lx = 60000; 
Ly = 13000;
ds = 50;
m = Lx/ds;
n = Ly/ds;
rolloff = 1/4000;
[z , PixelWidth, PSD] = artificial_surf(deltaD, H, Lx, m , n, rolloff); 

x = 0:ds:Lx-ds;
y = 0:ds:Ly-ds;
[X,Y] = meshgrid(x,y);
figure; mesh(X,Y,z)

rand_bed.X = X;
rand_bed.Y = Y;
rand_bed.z = z;

filename = ['dD',num2str(deltaD),'_H0',num2str(H*10),'_RO',num2str(1/rolloff),'.mat'];
%save(filename,'rand_bed')