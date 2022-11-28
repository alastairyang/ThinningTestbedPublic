% This script makes a scatter plot
% it quantifies 1. the magnitude of difference in dh/dt (MAE) 2.The differences
% in the shape (dynamic time warping)
ref_foldername = "analyzed_data/calve_only";
ref_folder_prefix = "ht_calve_";
foldernames = ["analyzed_data/mu_calve","analyzed_data/smw_calve","analyzed_data/gp3_calve"];
folder_prefixs = ["ht_mu_calve_","ht_smw_calve_","ht_gp3_calve_"];

symbs = [28,128,65; 229,87,9; 4,104,107]/255;

md_count = 0;
figure('Position',[100,100,600,500])

for i = 1:length(foldernames)
    ref_folder_dir = natsortfiles(dir([pwd '/' convertStringsToChars(ref_foldername)]));
    ref_folder_dir = struct2table(ref_folder_dir);
    folder_dir = natsortfiles(dir([pwd '/' convertStringsToChars(foldernames(i))]));
    folder_dir = struct2table(folder_dir);
    % remove  '.' and '..'
    bools = cellfun(@(s) ~strcmp(s(1),'.'), folder_dir.name);
    folder_dir = folder_dir(bools,:);
    bools = cellfun(@(s) ~strcmp(s(1),'.'), ref_folder_dir.name);
    ref_folder_dir = ref_folder_dir(bools,:);
    
    for j = 1:size(folder_dir,1)
        % skip the irrelevant ones
        md_count = md_count + 1;

        % reference model (calving-only)
        ref_md = load(string(ref_folder_dir.folder(j))+"/"+ string(ref_folder_dir.name(j))).ht_data;
        ref_modelname = ref_folder_dir.name(j);
        ref_modelname = ref_modelname{1}(length(convertStringsToChars(ref_folder_prefix))+1:end-4);

        % model of interest
        md = load(string(folder_dir.folder(j))+"/"+ string(folder_dir.name(j))).ht_data;
        modelname = folder_dir.name(j);
        modelname = modelname{1}(length(convertStringsToChars(folder_prefixs(i)))+1:end-4);

        % h(t) -> (h(t) - h0)
        md.h = md.h - md.h(1,:);
        ref_md.h = ref_md.h - ref_md.h(1,:);
        % find dh/dt
        smooth_window = 20;
        % to smooth data, for both the reference model and the target model
        for line_i = 1:size(md.h,2)
            h_smooth(:,line_i) = smooth(md.h(:,line_i), smooth_window);
        end
        dh = diff(h_smooth,1,1);
        dt = md.t(2) - md.t(1);
        md_dhdt = dh/dt;
        % ref_md
        for line_i = 1:size(ref_md.h,2)
            h_smooth(:,line_i) = smooth(ref_md.h(:,line_i), smooth_window);
        end
        dh = diff(h_smooth,1,1);
        ref_md_dhdt = dh/dt;


        % measure MAE and dtw
        MAEs = mean(abs(md_dhdt - ref_md_dhdt));
        % normalize timeseries, then measure the 
        md_h_norm = md_dhdt./max(abs(md_dhdt));
        ref_md_h_norm = ref_md_dhdt./max(abs(ref_md_dhdt));
        dtw_dists = [];
        for k = 1:size(md_h_norm,2)
%             coef = corrcoef(md_h_norm(:,k), ref_md_h_norm(:,k));
%             coef = coef(1,2);
            dtw_dists = [dtw_dists, dtw(md_h_norm(:,k), ref_md_h_norm(:,k))];
        end
%         MAE = mean(MAEs);
%         dtw_dist = max(dtw_dists);

        scatter(MAEs, dtw_dists, 25,symbs(i,:),'filled','MarkerFaceAlpha',0.4);
        hold on
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
xlabel('Mean Absolute Error','Interpreter','latex','FontSize',13)
ylabel('Shape dissimilarity','Interpreter','latex','FontSize',13)

saveas(gcf, 'plots/mu_smw_distance.pdf')
