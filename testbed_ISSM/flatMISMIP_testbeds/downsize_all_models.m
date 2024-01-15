%% Reduce the size of all the ISSM models to be uploaded to repo
% Author: Donglai Yang
% Date: July 3, 2023

%% Deep downsizing
% This means we not only remmove extra TransientSolution and friction
% coefficient array, we also remove columns in the TransientSolution we are keeping
% This applies to all models, except the effective pressure experiment
deep_flag = 1;
modelname_skip = ["C0.mat","H0.mat","Zb.mat","k_budd.mat","MISMIP_yangTransient_Calving_MassUnloading.mat","MISMIP_yangTransient_CalvingOnly.mat"];
% find all model directories 
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);

for folder_i = 16:size(foldernames_tbl,1)
    foldername = [foldernames_tbl.folder{folder_i} '/' foldernames_tbl.name{folder_i} '/*.mat'];
    modelnames = struct2table(dir(foldername));
    % iterate over model names which we want to remove from the downsizing
    % process
    keep = ones(size(modelnames,1),1);
    for jj = 1:length(modelname_skip)
        keep_j = ~strcmp(string(modelnames.name), modelname_skip(jj));
        keep = keep.*keep_j;
    end
    keep = logical(keep);

    modelnames = modelnames(keep,:);
    disp(['Model: ' foldernames_tbl.name{folder_i}])
    for md_i = 1:size(modelnames,1) % 1:size(modelnames,1)
        disp(['     Working on ' modelnames.name{md_i}])
        modeldir = [modelnames.folder{md_i} '/' modelnames.name{md_i}];
        load(modeldir)
        % in "downsize_md" function, see the fields we'd like to remove
        md = downsize_md(md, deep_flag);
        % replace the original model file
        save(modeldir, 'md')
        disp('     Downsized model saved!')

    end

end
%% shallow downsizing
% For effective pressure experiment, we do not remove any columns in the
% TransientSolution. We merely remove any extra TransientSolution fields
% and friction coefficient array (hence "shallow downsizing")

deep_flag = 0;
modelname_keep = 'MISMIP_yangTransient_Calving_MassUnloading.mat';
% find all model directories 
foldernames = natsortfiles(dir([pwd,'/long_models_yang']));
foldernames_tbl = struct2table(foldernames);
bools = cellfun(@(s) ~strcmp(s(1),'.'), foldernames_tbl.name);
foldernames_tbl = foldernames_tbl(bools,:);

% iterate over testbeds
for folder_i = 1:size(foldernames_tbl,1)
    foldername = [foldernames_tbl.folder{folder_i} '/' foldernames_tbl.name{folder_i} '/*.mat'];
    modelnames = struct2table(dir(foldername));
    keep_i = strcmp(string(modelnames.name), modelname_keep);
    modelnames = modelnames(keep_i,:);
    disp(['Model: ' foldernames_tbl.name{folder_i}])
    % iterate over experiments of a testbed
    for md_i = 1:size(modelnames,1)
        disp(['     Working on ' modelnames.name{md_i}])
        modeldir = [modelnames.folder{md_i} '/' modelnames.name{md_i}];
        load(modeldir)
        md = downsize_md(md, deep_flag);
        % replace the original model file
        save(modeldir, 'md')
        disp('     Downsized model saved!')

    end

end
