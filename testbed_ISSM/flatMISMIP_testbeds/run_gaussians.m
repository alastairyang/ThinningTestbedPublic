% This script runs all 4 gaussian perturbation experiments with different
% magnitudes of basal shear stress reduction
gauss_mag_vec = [0.1, 0.2, 0.3, 0.4]; % 10, 20, 30, 40 percent
simu_group = 1;

%% start

% iterate over the magnitudes
for k = 2:4
    
    % first, scan the current directory and delete all previously partitioned
    % run scripts
    i = 0;
    while true
        i = i + 1;
        filename = 'runme_gauss.m';
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
    fprintf(fid, ['function runme_gauss()',' \n']);
    fprintf(fid, ['md_idx = ', '[',num2str(current_idx),']','; \n']);
    fprintf(fid, ['tbl_filename = "runtime_table_gauss_', num2str(k),'.csv"','; \n']);
    fprintf(fid, ['gauss_mag = ', num2str(gauss_mag_vec(k)),'; \n']);
    fwrite(fid, frest);
    fprintf(fid, ['end','\n']);

    % run the model
    runme_gauss

    % delete filename
    delete('runme_gauss.m')
end
    
