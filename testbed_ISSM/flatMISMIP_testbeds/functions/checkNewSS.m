function [rel_GL, front, dH, t_plot, x, y, gls_expt] = checkNewSS(md_ctrl, md_expt, ds, nt)
% CheckNewSS: Check when new steady state is reached (thickness change approaching zero)
% after the perturbation stops.

    % parameter
    dGL_tol = 1; %  m/a

    % main
    expt = struct2table(md_expt.results.TransientSolution);
    ctrl = struct2table(md_ctrl.results.TransientSolution);
    elements = md_ctrl.mesh.elements;

    % relative grounding line
    gls_expt = zeros(size(expt.time));
    gls_ctrl = zeros(size(ctrl.time));
    front = zeros(size(ctrl.time));
    for i = 1:length(expt.time); gls_expt(i) = locate_groundingline(md_expt, md_expt.results.TransientSolution(i).MaskOceanLevelset); end
    for i = 1:length(ctrl.time); gls_ctrl(i) = locate_groundingline(md_ctrl, md_ctrl.results.TransientSolution(i).MaskOceanLevelset); end
    for i = 1:length(ctrl.time); front(i) = locate_calvingfront(md_ctrl, md_ctrl.results.TransientSolution(i).MaskIceLevelset); end

    % find which simulation is shorter
    if ctrl.time(end) < expt.time(end)
        expt_H_interp = transpose(interp1(expt.time, [expt.Thickness{:}]', ctrl.time));
        deltaH = expt_H_interp - [ctrl.Thickness{:}];
        deltaH_cell = num2cell(deltaH,1);
        [md_grid, ~, ~] = mesh_to_grid_overtime(elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
        % mask out non-ice part
        [mask_grid, x, y] = mesh_to_grid_overtime(elements, md_ctrl.mesh.x, md_ctrl.mesh.y, ctrl.MaskIceLevelset, ds);
        md_grid(mask_grid>0) = nan;
        dH = permute(md_grid,[2,3,1]);

        % relative grounding line
        gls_expt = interp1(expt.time, gls_expt, ctrl.time);
        rel_GL = gls_expt - gls_ctrl;

        t = ctrl.time - ctrl.time(1);

    else % expt has a shorter simulation time span
        ctrl_H_interp = transpose(interp1(ctrl.time, [ctrl.Thickness{:}]', expt.time));
        deltaH = [expt.Thickness{:}] - ctrl_H_interp;
        deltaH_cell = num2cell(deltaH,1);
        [md_grid, x, y] = mesh_to_grid_overtime(elements, md_ctrl.mesh.x, md_ctrl.mesh.y, deltaH_cell, ds);
        % mask out non-ice part
        [mask_grid, ~, ~] = mesh_to_grid_overtime(elements, md_ctrl.mesh.x, md_ctrl.mesh.y, expt.MaskIceLevelset, ds);
        md_grid(mask_grid>0) = nan;
        dH = permute(md_grid,[2,3,1]);

        % relative grounding line
        gls_ctrl = interp1(ctrl.time, gls_ctrl, expt.time);
        rel_GL = gls_expt - gls_ctrl;
        % front
        front = interp1(ctrl.time, front, expt.time);

        t = expt.time - expt.time(1);

    end

    % get center flow line data
    mid_i = floor(size(dH,1)/2);
    dH_mids = squeeze(dH(mid_i,:,:));
    % when the thickness change stops
%     max_H = max(abs(dH_mids), [], 1);
%     max_dHdt = transpose((max_H(2:end) - max_H(1:end-1)))./(t(2:end)-t(1:end-1));
    % when the relative GL stops moving
    %rel_GL_rate = (rel_GL(2:end)-rel_GL(1:end-1))./(t(2:end)-t(1:end-1));
    expt_GL_rate = (gls_expt(2:end) - gls_expt(1:end-1))./(t(2:end)-t(1:end-1));
%     idx_dH      = find(abs(max_dHdt)<dH_tol,1,'first');
    %idx_rel_GL  = find(abs(rel_GL_rate)<dGL_tol,1,'first');
    idx = find(abs(expt_GL_rate)<dGL_tol,1,'first');
    
%    idx = max([idx_dH, idx_rel_GL, idx_expt_GL]);
    t_plot = linspace(t(1),t(end),nt); % skip 20 years since the transientrestart can cause numerical error that acts effectively like a tiny forcing

    % interpolate dH and relative GL at plus minus 5 years around this time
    dH       = transpose(interp1(t, dH_mids', t_plot));
    rel_GL   = interp1(t, rel_GL, t_plot');
    gls_expt = interp1(t,gls_expt, t_plot');
    front    = interp1(t, front, t_plot');


%Previous solution: put a n-year window around the steady state time
%     %    idx = max([idx_dH, idx_rel_GL, idx_expt_GL]);
%     if isempty(idx) 
%         warning('No steady state reached!');
%         t_stop = t(end)-50;
%     else
%         t_stop = t(idx); 
%     end
%     t_plot = t_stop:dt:t_stop+t_window;
% 
%     % interpolate dH and relative GL at plus minus 5 years around this time
%     dH       = transpose(interp1(t, dH_mids', t_plot));
%     rel_GL   = interp1(t, rel_GL, t_plot');
%     gls_expt = interp1(t,gls_expt, t_plot');
%     front    = interp1(t, front, t_plot');


    
end
