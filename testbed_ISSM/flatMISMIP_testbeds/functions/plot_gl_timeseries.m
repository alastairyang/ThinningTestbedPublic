function [gl_locs,t] = plot_gl_timeseries(md, make_plot)
%PLOT_GL_TIMESERIES Summary of this function goes here

    nt = size(md.results.TransientSolution,2);
    gl_locs = zeros(1,nt);
    t = zeros(1,nt);
    for k = 1:nt
        gl_mask = md.results.TransientSolution(k).MaskOceanLevelset;
        gl_locs(k) = locate_groundingline(md, gl_mask);
        t(k) = md.results.TransientSolution(k).time;
    end

    if nargin > 1 && make_plot
        figure; plot(1:nt, gl_locs,'-r','LineWidth',2);
    end
    fprintf('Grounding line location acquired!\n')

end

