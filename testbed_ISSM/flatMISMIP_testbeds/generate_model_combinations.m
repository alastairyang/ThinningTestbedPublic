% Generate all combinations of model variable values 
% Current (June 19, 2022) variables of interest include:
%   1. Fjord width
%   2. Grounding line depth
%   3. Basal sliding
%   4. Background basal friction level

%% Specify low level and high level
% [low, high]
fjord_width = [5000, 8000, 11000]; % meter
gl_depth = [-400, -200, -100]; % meter, grouning line depth
bs_law = [1]; % p number in the Paterson sliding law, 1 is linear, 3 is nonlinear
bg_friccoef = [20, 40, 80]; % (unit?), background basal friction level

var_names = ["fjord_width","delta_groundingline_depth","basalfric_law","background_friccoef"];

%% Generate all combinations and export to a .csv spreadsheet
N = length(fjord_width)*length(gl_depth)*length(bs_law)*length(bg_friccoef);
model_vars = table('Size', [N, length(var_names)],...
                   'VariableTypes', {'single','single','int32','single'},...
                   'VariableNames', var_names);

% combinations
[A, B, C, D]=ndgrid(fjord_width, gl_depth, bs_law, bg_friccoef);
A1=reshape(A, [], 1);
B1=reshape(B, [], 1);
C1=reshape(C, [], 1);
D1 = reshape(D, [], 1);
d = [A1,B1,C1, D1];

% add to the table
model_vars{:,:} = d;
writetable(model_vars, 'md_var_combinations.csv')