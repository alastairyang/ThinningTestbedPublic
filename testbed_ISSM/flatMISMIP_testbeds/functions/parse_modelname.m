function [width, gl, fc] = parse_modelname(name)
%PARSE_MODELNAME This function returns the width, grounding line depth, and
%basal drag level from the model name alone

    % first split by '_'
    parts = split(name,'_');
    if strcmp(parts{end}(end-3:end), '.mat')
        % this modelname contains the extension; remove
        parts{end} = parts{end}(1:end-4);
    end
    % it's always the last three cells
    % 1. width, 2.grounding line depth, 3.basal friction level
    width = str2double(parts{end-2}(2:end));
    gl = str2double(parts{end-1}(3:end));
    fc = str2double(parts{end}(3:end));
    
end

