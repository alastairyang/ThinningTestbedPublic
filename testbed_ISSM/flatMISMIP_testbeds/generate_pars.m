%% This file generates .par files (to parameterize the models)
% The numeric values of certain geometric variables (e.g., fjord width) are
% from the "md_var_combination" file.

% This file generates .par scripts, which can be run in ISSM to create the
% geometry and parameterize the model.

%% Parameter
sens_run = 0; % no sensitivity run

%% Read in model combination data and create labels
mdvar_combs = readtable('md_var_combinations.csv');
N_var = size(mdvar_combs,2);
N_mds = size(mdvar_combs,1);
varnames = mdvar_combs.Properties.VariableNames;
vartypes = repelem(["int32"], length(varnames)); % repeat element
% create a new table for labels
mdvar_combs_label = table('Size', size(mdvar_combs),...
                          'VariableTypes', vartypes,...
                          'VariableNames', varnames);

% create dictionary for the variables and their label initials
initial_dict = create_initial_dict();
% get the sequence of the label initials
varname_labels = cellfun(@(key) initial_dict(key), varnames);

% creating labels
% lower value is 0, high values: just keep adding 1 to the previous one
% e.g., if only two values, [0,1]; three values, [0,1,2]
for i = 1:length(varnames)
    col_data = mdvar_combs.(varnames{i});
    col_data_uniq = unique(col_data);
    M = containers.Map(transpose(sort(col_data_uniq)), 0:(length(col_data_uniq)-1));
    col_data_idx = arrayfun(@(key) M(key), col_data);
    mdvar_combs_label.(varnames{i}) = col_data_idx;    
end

% combine index and initials to form labels
% these few lines are so elegant. I won't be able to come up on my own
mdvar_combs_label = num2str(mdvar_combs_label.Variables,'%d');
varname_labels_matrix = repmat(varname_labels, [N_mds,1]);
varname_labels_matrix = varname_labels_matrix';
mdvar_combs_label = mdvar_combs_label';
row_interleave = reshape([varname_labels_matrix(:) mdvar_combs_label(:)]',2*size(varname_labels_matrix,1), []);
labels = row_interleave(:)';
% now labels is a long string, we need to use regular expression to split
labels = regexp(labels, sprintf('\\w{1,%d}', N_var*2), 'match');

%% Generate the .par files
for i = 1:N_mds
    label = labels(i);
    label = label{:};
    var_table = mdvar_combs(i,:);
    var_table.('label') = label;
    write_par(var_table)
end

%% FUNCTION:
function dict = create_initial_dict()
    keys = {'fjord_width', 'delta_groundingline_depth', 'basalfric_law', 'background_friccoef'};
    vals = {'f',            'g',                 'l',             'b'};
    dict = containers.Map(keys, vals);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function to generate .par files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_par(var_table)

    Ly = var_table.('fjord_width');
    delta_gl_depth = var_table.('delta_groundingline_depth');
    bg_fric_coef = var_table.('background_friccoef');

    par_name = ['par','_W', num2str(Ly), '_GL', num2str(delta_gl_depth), '_FC', num2str(bg_fric_coef)];
    fid = fopen('par_files/MISMIP_template.par');
    % remove the first four lines
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    % read in the rest
    frest = fread(fid, Inf);
    fclose(fid);
    % destination file
    fid = fopen(['par_files/', par_name,'.par'], 'w');
    fprintf(fid, ['Ly = ', num2str(Ly),'; \n']);
    fprintf(fid, ['delta_gl_depth = ', num2str(delta_gl_depth), '; \n']);
    fprintf(fid, ['bg_fric_coef = ', num2str(bg_fric_coef),'; \n']);
    fwrite(fid, frest);

end