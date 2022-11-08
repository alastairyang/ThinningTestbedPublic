% this script creates multiple versions runme_free_front_long_*.m, with
% mutually exclusive independent model runs assigned. e.g., runme_..._part_1, ...part_2, part_3,...

% first, scan the current directory and delete all previously partitioned
% run scripts
i = 0;
while true
    i = i + 1;
    filename = ['runme_part_',num2str(i),'.m'];
    if exist(filename,'file')
        disp(['delete ' filename])
        delete(filename);
    else
        disp([filename ' not found! Break out.'])
        break;
    end
end

% create partitions
N_files = input('How many partitions?\n');
model_idx = cell(1,N_files);
% request model indices for each partitioned run
for j = 1:N_files
    prompt = ['model indices for ', num2str(j),'\n'];
    current_idx = input(prompt);

    % check if these indices overlap with previous ones
    previous_model_idx = cell2mat(model_idx);
    if ~isempty(intersect(previous_model_idx, current_idx))
        error('Model indices are not mutually exclusive! Repeated index found')

    else
        model_idx{j} = current_idx;
        % duplicate runme file for partitioned run
        fid = fopen('runme_free_front_long.m','r');
        % skip the first n lines, n = ? Check the base template for details
        n = 6;
        for i = 1:n
            fgetl(fid); 
        end
        % read in the rest
        frest = fread(fid, Inf);
        fclose(fid);
        % destination file
        fid = fopen(['runme_part_', num2str(j),'.m'], 'w');
        fprintf(fid, ['md_idx = ', '[',num2str(model_idx{j}),']','; \n']);
        fprintf(fid, ['tbl_filename = "runtime_table_', num2str(j),'.csv"','; \n']);
        fwrite(fid, frest);
    end
end

