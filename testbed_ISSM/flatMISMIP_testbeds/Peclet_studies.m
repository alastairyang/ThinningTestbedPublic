% peclet number studies
nt = size(md.results.TransientSolution,2);
ds = 50;
smooth_span = 30;
for i = 1
    year = md.results.TransientSolution(i).time;
    index = md.mesh.elements;
    meshx = md.mesh.x;
    meshy = md.mesh.y;
    k = mesh_to_grid(index, meshx, meshy, md.friction.C(1:end-1,i));
    H = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Thickness);
    hw = md.results.TransientSolution(i).Base;
    hw(hw>0) = 0;
    pw = hw.*(md.materials.rho_water/md.materials.rho_ice);
    pw = mesh_to_grid(index, meshx, meshy, pw);
    s = mesh_to_grid(index, meshx, meshy, md.results.TransientSolution(i).Surface);
    alpha = zeros(size(s));
    for j = 1:size(s,1)
        alpha(j,:) = smooth(gradient(s(j,:), ds), 30);
    end
    % now calculate C0, D0, grad(D0)
    C0 = compute_C0(k, H, alpha, pw);
    D0 = compute_D0(k, H, alpha, pw);
    grad_D0 = zeros(size(D0));
    for j = 1:size(s,1)
        grad_D0(j,:) = smooth(gradient(D0(j,:), ds), 30);
    end
    %grad_D0 = compute_grad_D0(D0, ds);
    % distance to ice front
    distance = mesh_to_grid(index, meshx, meshy, -1*md.results.TransientSolution(i).MaskIceLevelset);
    distance(distance<0) = 0;
    Pe = (C0 - grad_D0)./D0.*distance;
    % Pe across-flow average
    % the middle 20
    if rem(size(Pe,1), 2) == 0
        mid_i = size(Pe,1)/2;
    else
        mid_i = (size(Pe,1)+1)/2;
    end
    Pe_mean = mean(Pe(mid_i-10:mid_i+10,:),1,'omitnan');
    
end
%% plot
x = 0:ds:(size(Pe_mean,2)-1)*ds;
figure;
plot(x, low_fric_Pe);hold on
plot(x, high_fric_Pe); hold off
ylim([0,5])
xlim([10000,50000])
legend(["low friction","high friction"])
xlabel('along flow distance (m)')
ylabel('Peclet number')
exportgraphics(gcf,'Pe_lowhigh_fric.png','Resolution',300)

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