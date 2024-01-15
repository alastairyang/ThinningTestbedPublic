function bool = compare_GLvalue(str, val)
%   This function extracts the grounding line value from a model run
%   filename and compares it to a supplied numerical value. It outputs a
%   boolean value: if identical, it outputs 1
    filename_split = split(str, '_');
    initials = string(cellfun(@(s) s(1:2), filename_split, 'UniformOutput', false));
    GLvalue = filename_split(strcmp('GL', initials));
    GLvalue = str2double(GLvalue{1}(3:end));

    bool = GLvalue == val;
end