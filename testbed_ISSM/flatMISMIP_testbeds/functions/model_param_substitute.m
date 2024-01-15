function md = model_param_substitute(md, sens_data, x, y)
%MODEL_PARAM_SUBSTITUTE This function substitutes parameter fields of the ISSM model
% This can include forcing data or other parameters
%
% Input:
%       md[model]      : ISSM model
%       sub_data[array]: cell array (for time-dependent data) storing the data for substitution
%       x[array]       : a double vector containing x axis
%       y[array]       : a double vector containing y axis
%
% Output:
%       md[model]      : return the modified ISSM model

    data = sens_data.data;
    nt = numel(data);
    fieldname = sens_data.var;

    if strcmp('shelf_melt', fieldname)
        all_melt_rates = [];
        % convert the cell array into one matrix
        for i = 1:nt
            all_melt_rates = [all_melt_rates, data{i}];
        end
        md.basalforcings.floatingice_melting_rate = InterpFromGridToMesh(x',y,all_melt_rates,md.mesh.x,md.mesh.y,0);
    end

    if strcmp('fric_coef', fieldname)
        all_fric_coef = [];
        % convert the cell array into one matrix
        for i = 1:nt
            all_fric_coef = [all_fric_coef, data{i}];
        end
        md.friction.coefficient = InterpFromGridToMesh(x',y,all_fric_coef,md.mesh.x,md.mesh.y,0);
    end

    if strcmp('rheoB', fieldname)
        all_rheoB = [];
        % convert the cell array into one matrix
        for i = 1:nt
            all_rheoB = [all_rheoB, data{i}];
        end
        md.materials.rheology_B = InterpFromGridToMesh(x',y,all_rheoB,md.mesh.x,md.mesh.y,0);
    end
end
