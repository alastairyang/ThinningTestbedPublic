%% Ice volume timeseries
% we are interested in the comparison between deep - shallow grounding
% line, with - without the effective pressure feedback
%% parameters
foldername = 'long_models_yang';

calving_modelname = 'MISMIP_yangTransient_CalvingOnly.mat';
calving_mu_modelname = 'MISMIP_yangTransient_Calving_MassUnloading.mat';

shal_name = 'model_W11000_GL0_FC120000';
deep_name = 'model_W11000_GL400_FC120000';

%% load
% NOTE!! If you change the model paths below, you also need to change the
% legend strings and the associated color and style arrays
md1_path = string([foldername,'/', shal_name,'/',calving_modelname]);
md2_path = string([foldername,'/', shal_name,'/',calving_mu_modelname]);
md3_path = string([foldername,'/', deep_name,'/',calving_modelname]);
md4_path = string([foldername,'/', deep_name,'/',calving_mu_modelname]);
md_paths = [md1_path, md2_path, md3_path, md4_path];
legend_strs = ["shallow, retreat only","shallow, retreat + feedback",...
               "deep, retreat only",   "deep, retreat + feedback"];
colors = [252,175,124;
          252,175,125;
          135,201,195;
          135,201,195]/255;
styles = ["-","-","-","-"];
markers = ["none","^","none","^"];

%% extract data
times = cell(1,length(md_paths));
ice_volumes_grad = cell(1,length(md_paths));
ice_volumes = cell(1,length(md_paths));
for i = 1:length(md_paths)
    md = load(md_paths(i)).md;
    results_tbl = struct2table(md.results.TransientSolution);
    times{i} = transpose(results_tbl.time);
    dt = mean(times{i}(2:end) - times{i}(1:end-1));
    ice_volumes_grad{i} = transpose(gradient(smooth(results_tbl.IceVolume, 20),0.1));
    ice_volumes{i} = transpose((results_tbl.IceVolume));
end

%% plot
figure('Position',[100,100,600,600])
p1 = subplot(2,1,1);
% ice volume
for i = 1:length(md_paths)
    plot(times{i}(1:10:end), ice_volumes{i}(1:10:end),...
         "Color", colors(i,:), 'LineStyle', styles(i),'Marker',markers(i),'LineWidth',2);
    hold on
end
set(gca,'XTickLabel',[]);
% create path to fill in the area between curves with and without feedback
patch([times{1} fliplr(times{2})], [ice_volumes{1} fliplr(ice_volumes{2})], colors(1,:))
patch([times{3} fliplr(times{4})], [ice_volumes{3} fliplr(ice_volumes{4})], colors(3,:))
% the retreat stops at the 22nd year in the perturbation. We grey out the
% years beyond
yr_stop = times{1}(1) + 22;
[d,ix] = min(abs(times{1}-yr_stop));
area(times{4}(ix:end),5.5e11*ones(size(times{4}(ix:end))),0)
alpha(0.2)
hold off;
ylim([3.4e11,4.8e11])
xlim([times{4}(1), times{4}(end)])
ylabel('Ice volume ($m^3$)','Interpreter','latex','FontSize',16)

% time rate of ice volume
p2 = subplot(2,1,2);
for i = 1:length(md_paths)
    plot(times{i}(1:10:end), ice_volumes_grad{i}(1:10:end),...
         "Color", colors(i,:), 'LineStyle', styles(i),'Marker',markers(i),'LineWidth',2);
    hold on
end
% create path to fill in the area between curves with and without feedback
patch([times{1} fliplr(times{2})], [ice_volumes_grad{1} fliplr(ice_volumes_grad{2})], colors(1,:))
patch([times{3} fliplr(times{4})], [ice_volumes_grad{3} fliplr(ice_volumes_grad{4})], colors(3,:))
% the retreat stops at the 22nd year in the perturbation. We grey out the
% years beyond
yr_stop = times{1}(1) + 22;
[d,ix] = min(abs(times{1}-yr_stop));
area(times{4}(ix:end),-9e9*ones(size(times{4}(ix:end))),2e9)
alpha(0.2)
hold off;
ylim([-8e9,0])
xlim([times{4}(1), times{4}(end)])
q1 = get(p1, 'Position');
q2 = get(p2, 'Position');
q1(2) = q2(2)+q2(4);
set(p1, 'pos', q1);
xlabel('Time (yr)','Interpreter','latex','FontSize',16)
ylabel('Ice volume change rate ($m^3/a$)','Interpreter','latex','FontSize',16)
legend(legend_strs,'Location','southwest','box','off')

saveas(gcf, 'plots/ice_volume.pdf')
