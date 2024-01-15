function [s_profiles, base_profiles, bed_profile, t] = sample_profile_evol(md, dt, ds, t_interval)
%SAMPLE_PROFILE_EVOL sample the lateral profile of glacier (along the
%center flow line) at multiple timesteps

    % define sampled time slice indices
    nt = size(md.results.TransientSolution,2);
    skip_i = floor(t_interval/dt);
    ti_sampled = 1:skip_i:nt;

    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)
    
    % get along flow axis in regular meshgrid
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    x = 0:ds:Lx-ds;
    y = 0:ds:Ly-ds;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end

    % get bed profile which does not change overtime
    bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.geometry.bed,...
            x, y, NaN);
    bed_profile = bed(mid_i,:);

    % pre-allocate
    s_profiles = zeros(length(ti_sampled), length(bed_profile));
    base_profiles = zeros(length(ti_sampled), length(bed_profile));
    % iterate
    for i = 1:length(ti_sampled)
        ti = ti_sampled(i);
        surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(ti).Surface,...
            x, y, NaN);
        base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
            md.results.TransientSolution(ti).Base,...
            x, y, NaN);
        s_profiles(i,:)  = surface(mid_i,:);
        base_profiles(i,:) = base(mid_i,:);
        
    end
    % construct time
    t = ti_sampled*dt;
end

