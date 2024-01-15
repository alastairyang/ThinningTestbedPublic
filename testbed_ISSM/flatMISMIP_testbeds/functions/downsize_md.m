function md = downsize_md(md, deep_flag)
%DOWNSIZE_MD reduces the ISSM model file size by removing any extra data
%array
%
%   Input:
%       md [ISSM model]: ISSM model class
%       deep_flag [bool]: 1 - deep cleaning; 0 or empty: no deep cleaning
%                        'deep cleaning' means removing the redundant columns
%                        in md.results.TransientSolution
%   
%   Output:
%       md [ISSM model]: downsized ISSM model

    % first show the original size
    dt = whos('md');
    dt_mb = dt.bytes*9.53674e-7;
    disp(['The orignal size is ' num2str(dt_mb) ' MB'])

    % downsize results
    if ~isempty(fieldnames(md.results)) % it has transient solution
        solution = md.results.TransientSolution;
        md.results = [];
        md.results.TransientSolution = solution;
    end
    
    % downsize friction array
    if size(md.friction.C, 2)>1 
        % we keep the last column
        md.friction.C = md.friction.C(:,end);
    end

    if deep_flag && ~isempty(fieldnames(md.results))
        keys = ["StrainRatexx",...
                "StrainRateyy",...
                "StrainRatexy",...
                "StrainRateeffective",...
                "DeviatoricStressxx",...
                "DeviatoricStressxy",...
                "DeviatoricStressyy",...
                "Pressure",...
                "Vx",...
                "Vy"];
        results_tbl = struct2table(md.results.TransientSolution);
        % iteratively remove variable columns
        for key = keys
            try
                results_tbl = removevars(results_tbl, key);
            catch
                disp(key + " not found! We skip it...")
            end
        end
        % assign back to ISSM md
        results = table2struct(results_tbl);
        md.results = []; md.results.TransientSolution = results;
    end

    % show the new size
    dt = whos('md');
    dt_mb = dt.bytes*9.53674e-7;
    disp(['The size after downsizing is ' num2str(dt_mb) ' MB'])
end

