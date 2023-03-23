load('long_models_yang/model_W11000_GL400_FC30000/MISMIP_yangTransient_Calving_MassUnloading.mat')

% peclet number studies
nt = size(md.results.TransientSolution,2);
ds = 50;
smooth_span = 30;
sampled_ti = 1:20:nt;
[~,x,y] = mesh_to_grid(index, meshx, meshy, md.friction.C(1:end-1,i),ds);
% pre-allocate
Pes = zeros(length(sampled_ti), length(x));

count = 0;
for i = sampled_ti
    count = count+1;
    year = md.results.TransientSolution(i).time;
    index = md.mesh.elements;
    meshx = md.mesh.x;
    meshy = md.mesh.y;
    [k,~,~] = mesh_to_grid(index, meshx, meshy, md.friction.C(1:end-1,i),ds);
    [H,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Thickness,ds);
    [mask,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).MaskOceanLevelset,ds);
    % make floating parts nan
    hw = md.results.TransientSolution(i).Base;
    hw(hw>0) = 0;
    pw = hw.*(md.materials.rho_water/md.materials.rho_ice);
    [pw,~,~] = mesh_to_grid(index, meshx, meshy, pw, ds);
    [s,~,~] = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Surface, ds);
    alpha = zeros(size(s));
    smooth_n = 2000/ds;
    for j = 1:size(s,1)
        alpha(j,:) = smooth(gradient(s(j,:), ds), smooth_n);
    end
    % now calculate C0, D0, grad(D0)
    C0 = compute_C0(k, H, alpha, pw);
    D0 = compute_D0(k, H, alpha, pw);
    grad_D0 = zeros(size(D0));
    for j = 1:size(s,1)
        grad_D0(j,:) = smooth(gradient(D0(j,:), ds), smooth_n);
    end
    %grad_D0 = compute_grad_D0(D0, ds);
    % distance to ice front
    [distance,~,~] = mesh_to_grid(index, meshx, meshy, -1*md.results.TransientSolution(i).MaskIceLevelset, ds);
    distance(distance<0) = 0;
    Pe = (C0 - grad_D0)./D0.*distance;
    % Pe across-flow average
    % the middle 20
    if rem(size(Pe,2), 1) == 0
        mid_i = size(Pe,1)/2;
    else
        mid_i = (size(Pe,1)+1)/2;
    end
    % only look at the middle 20 lines
    Pe(mask<0) = nan;
    Pes(count,:) = mean(Pe(mid_i-10:mid_i+10,:),1,'omitnan');
    
end
%% plot
plot_x = 0:ds:(size(Pes,2)-1)*ds;
figure; plot(plot_x/1000, Pes,'LineWidth',1);ylim([0,9]);colororder(cool(13))
xlabel('Along flow distance (km)','FontName','Aria','FontSize',15)
ylabel('Peclet Number','FontName','Aria','FontSize',15)
plotname = ['Pe_evolv_', md.miscellaneous.name(9:end),'.png'];
exportgraphics(gcf, ['plots/', plotname], 'Resolution',400)

%% functions
function C0 = compute_C0(k, H, alpha, pw)
    m = 1;
    n = 3;
    C0 = k.*((H.*alpha).^n./(H - pw).^m).*(1 + (n-m)*(H - n*pw)./(H - pw));
end

function D0 = compute_D0(k, H, alpha, pw)
    m = 1;
    n = 3;
    D0 = k.*n.*(alpha.^(n-1).*H.^(n+1)./((H - pw).^m));
end