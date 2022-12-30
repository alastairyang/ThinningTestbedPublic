% This script runs all 4 gaussian perturbation experiments with different
% magnitudes of basal shear stress reduction
gauss_mags = [0.1, 0.2, 0.3, 0.4]; % 10, 20, 30, 40 percent
simu_group = 1:18;

%% start
% this script creates multiple versions runme_free_front_long_*.m, with
% mutually exclusive independent model runs assigned. e.g., runme_..._part_1, ...part_2, part_3,...

% iterate over the magnitudes
for k = 1:length(gauss_mags)
    
    % first, scan the current directory and delete all previously partitioned
    % run scripts
    i = 0;
    while true
        i = i + 1;
        filename = ['runme_part_gauss_',num2str(i),'.m'];
        if exist(filename,'file')
            disp(['delete ' filename])
            delete(filename);
        else
            disp([filename ' not found! Break out.'])
            break;
        end
    end

    % create a run file 
    current_idx = simu_group;

    % duplicate runme file for partitioned run
    fid = fopen('runme_free_front_long.m','r');
    % skip the first n lines, n = ? Check the base template for details
    n = 7;
    for i = 1:n
        fgetl(fid);
    end
    % read in the rest
    frest = fread(fid, Inf);
    fclose(fid);
    % destination file
    fid = fopen('runme_gauss.m', 'w');
    fprintf(fid, ['md_idx = ', '[',num2str(current_idx),']','; \n']);
    fprintf(fid, ['tbl_filename = "runtime_table_', num2str(1),'.csv"','; \n']);
    fprintf(fid, ['gauss_mag = ', num2str(gauss_mags(k)),'; \n']);
    fwrite(fid, frest);

    % run the model
    runme_gauss

    % delete filename
    delete('runme_gauss.m')
end
    
