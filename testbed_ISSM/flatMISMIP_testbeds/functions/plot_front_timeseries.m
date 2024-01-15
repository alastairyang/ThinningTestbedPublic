function [front_locs,t] = plot_front_timeseries(md, make_plot)
%PLOT_FRONT_TIMESERIES get the location of the front over time

    nt = size(md.results.TransientSolution,2);
    t = zeros(1,nt);
    front_locs = zeros(1,nt);
    for k = 1:nt
        front_mask = md.results.TransientSolution(k).MaskIceLevelset;
        front_locs(k) = locate_calvingfront(md, front_mask);
        t(k) = md.results.TransientSolution(k).time;
    end
    
    if nargin > 1 && make_plot
        figure; plot(1:nt, front_locs,'-r','LineWidth',2);
    end
    fprintf('Frontal location acquired!\n')
end

