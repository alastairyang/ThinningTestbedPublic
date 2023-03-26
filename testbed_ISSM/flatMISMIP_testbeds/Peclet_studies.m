load('long_models_yang/model_W5000_GL400_FC30000/MISMIP_yangTransient_Calving_MassUnloading.mat')

% peclet number studies
index = md.mesh.elements;
meshx = md.mesh.x; meshy = md.mesh.y;
nt = size(md.results.TransientSolution,2);
ds = 50;
smooth_span = 30;
dt = 0.1; % dt in the simulation results
t_interval = 2; % sample every 2 year
ti_interval = floor(t_interval/dt);
sampled_ti = 1:ti_interval:nt;
[~,x,y] = mesh_to_grid(index, meshx, meshy, md.friction.C(1:end-1,1),ds);

% pre-allocate
Pes = zeros(length(sampled_ti), length(x));
D0s = zeros(length(sampled_ti), length(x));
dHs = zeros(length(sampled_ti), length(x));

% thickness and basal coefficient at the steady state
[H0,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(1).Thickness,ds);
[k0,~,~] = mesh_to_grid(index, meshx, meshy, md.friction.C(1:end-1,1),ds);

count = 0;
for i = sampled_ti
    count = count+1;
    year = md.results.TransientSolution(i).time;
    index = md.mesh.elements;
    meshx = md.mesh.x;
    meshy = md.mesh.y;
    [H,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Thickness,ds);
    [mask,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).MaskOceanLevelset,ds);
    % make floating parts nan
    hw = md.results.TransientSolution(i).Base;
    hw(hw>0) = 0;
    pw = hw.*(md.materials.rho_water/md.materials.rho_ice);
    [pw,~,~] = mesh_to_grid(index, meshx, meshy, pw, ds);
    [s,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Surface, ds);
    alpha = zeros(size(s));
    smooth_n = 4500/ds; % Dennis used 10 times the ice thickness (here assuming h = 450)
    for j = 1:size(s,1)
        alpha(j,:) = smooth(gradient(s(j,:), ds), smooth_n);
    end
    % now calculate C0, D0, grad(D0)
    C0 = compute_C0(k0, H, alpha, pw);
    D0 = compute_D0(k0, H, alpha, pw);
    grad_D0 = zeros(size(D0));
    for j = 1:size(s,1)
        grad_D0(j,:) = smooth(gradient(D0(j,:), ds), smooth_n);
    end
    % distance to ice front is the characteristic length in Peclet number
    % calculation
    [distance,~,~] = mesh_to_grid(index, meshx, meshy, -1*md.results.TransientSolution(i).MaskIceLevelset, ds);
    distance(distance<0) = 0;
    Pe = (C0 - grad_D0)./D0.*distance;
    % Pe across-flow average
    % the middle 20 flowlines
    if rem(size(Pe,2), 1) == 0
        mid_i = size(Pe,1)/2;
    else
        mid_i = (size(Pe,1)+1)/2;
    end
    % only look at the middle 20 lines and behind the grounding line
    Pe(mask<0) = nan;
    Pes(count,:) = mean(Pe(mid_i-10:mid_i+10,:),1,'omitnan');
    D0s(count,:) = mean(D0(mid_i-10:mid_i+10,:),1,'omitnan');

    % get the thickness change from the steady state 
    dH = H-H0; dH(mask<0) = nan;
    dHs(count,:) = dH(mid_i,:);

end
%% Recover the glacier profiles overtime
[surface, base, bed, profiles_t] = sample_profile_evol(md, dt, ds, t_interval);

%% plot
% load colormap
cp = load('plots/colormap/lajolla.mat').lajolla;

% make the plot
plot_x = 0:ds:(size(Pes,2)-1)*ds;
figure('Position',[100,100,900,600])
tiledlayout(2,1,"TileSpacing","none")
nexttile
plot(x/1000, surface); hold on;
plot(x/1000, base); hold on;
plot(x/1000, bed,'-k'); hold off
ylim([-600,1300])
line_colors = colormap_to_colororder(cp, size(surface,1),1,100);
colororder(line_colors)

nexttile
yyaxis left % plot the Peclet number
plot(plot_x/1000, Pes,'LineWidth',1);
ylabel('Peclet Number','FontName','Aria','FontSize',15)
ylim([0,10]); set(gca,'XColor',[0 0 0]); set(gca,'YColor',[0 0 0]);
colororder(line_colors)
yyaxis right % plot dH
plot(plot_x/1000, dHs, '-.','LineWidth',1);
ylabel('dH (m)','FontName','Aria','FontSize',15)
xlabel('Along flow distance (km)','FontName','Aria','FontSize',15)
set(gca,'XColor',[0 0 0]); set(gca,'YColor',[0 0 0]);
colororder(line_colors)
plotname = ['Pe_evolv_', md.miscellaneous.name(9:end),'.png'];
exportgraphics(gcf, ['plots/', plotname], 'Resolution',400)

%% functions
function C0 = compute_C0(k, H, alpha, pw)
    m = 1;
    n = 3;
    C0 = k.*((H.*alpha).^n./(H - pw).^m).*(1 + ((n-m)*H - n*pw)./(H - pw));
end

function D0 = compute_D0(k, H, alpha, pw)
    m = 1;
    n = 3;
    D0 = k.*n.*(alpha.^(n-1).*H.^(n+1)./((H - pw).^m));
end