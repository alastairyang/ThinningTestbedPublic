%% parameters
foldername = 'long_models_yang';
model_type = 'MISMIP_yangTransient_Remeshexperiment.mat';

%% Visualize the Vx at the last time step
folder_dir = dir([pwd '/' foldername]);

for i = 1:size(folder_dir,1)
    % skip the irrelevant ones
    if ~strcmp(folder_dir(i).name(1), 'm')
        continue
    else
        md = load([folder_dir(i).folder '/' folder_dir(i).name '/' model_type]);
        md = md.md;
        result_num = size(md.results.TransientSolution, 2);
        model_name = folder_dir(i).name(7:end);
        disp(['The model is ' model_name])
        % print result
        plotmodel(md, 'data', md.results.TransientSolution(result_num).Vel,...
            'title',model_name);
        pause(1)

    end

end

%% Visualize an animated map of one specific model
% Visualize the Vx at the last time step
folder_dir = dir([pwd '/' foldername]);
i = 11;

md = load([folder_dir(i).folder '/' folder_dir(i).name '/' model_type]).md;
result_num = size(md.results.TransientSolution, 2);
model_name = folder_dir(i).name(7:end);
disp(['The model is ' model_name])
% print result
for j = 1:result_num
    plotmodel(md, 'data', md.results.TransientSolution(j).Vel,...
        'title',model_name, 'contourlevels',{0});
    %pause(0.2)
end



%% Visualize mean dh/dt timeseries
folder_dir = dir([pwd '/' foldername]);

figure;
plot_idx = 0;
for i = 1:size(folder_dir,1)
    % skip the irrelevant ones
    if ~strcmp(folder_dir(i).name(1), 'm')
        continue
    else
        plot_idx = plot_idx + 1;
        subplot(5,6, plot_idx)
        try
            md = load([folder_dir(i).folder '/' folder_dir(i).name '/' model_type]);
        catch
            continue
        end
        md = md.md;
        model_name = folder_dir(i).name(7:end);
        % calculating dh/dt
        try
            results_tbl = struct2table(md.results.TransientSolution);
        catch
            continue
        end
        result_num = size(md.results.TransientSolution, 2);
        dhdt = [];
        % calculate mean dh/dt
        for j = 2:result_num
            dhdt = [dhdt, mean(cell2mat(results_tbl.Thickness(j)) - cell2mat(results_tbl.Thickness(j-1)))/(results_tbl.time(j) - results_tbl.time(j-1))];
        end
        time = results_tbl.time;
        plot(time(1:end-1), dhdt)
        title(model_name)
        %ylim([-1,0])

    end

end

%% visualize dh/dt at various points on the center flowline
folder_dir = dir([pwd '/' foldername]);

figure;
plot_idx = 0;
md_count = 0;
for i = 1:size(folder_dir,1)
    % skip the irrelevant ones
    if ~strcmp(folder_dir(i).name(1), 'm')
        continue
    else
        md_count = md_count + 1;
        md = load([folder_dir(i).folder '/' folder_dir(i).name '/' model_type]).md;
        nt = size(md.results.TransientSolution,2);
        
        % geometry parameters
        Lx = max(md.mesh.x);
        Ly = max(md.mesh.y);
        ds = 250; % spacing, 250 meter
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        
        % the cell index of locations at the centerline
        sample_i = 35:5:80;
        thalweg_sample_ht = [];

        for j = 1:nt
            surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                                md.results.TransientSolution(j).Thickness,...
                                                x, y, NaN);
            surface_profile  = surface(mid_i,length(x) - sample_i);
            thalweg_sample_ht = [thalweg_sample_ht; surface_profile];
        end
        
         % time vector
         results_tbl = struct2table(md.results.TransientSolution);
         time = results_tbl.time;
         time_disp = [time(2:end);nan];
         dt = time_disp - time;
         dt = dt(1:end-1);
         % calculate dh/dt
         thalweg_sample_ht_disp = [thalweg_sample_ht(2:end,:);nan*ones(1,numel(sample_i))];
         dh = thalweg_sample_ht_disp - thalweg_sample_ht;
         dh = dh(1:end-1,:);
         dt_mtx = repmat(dt, 1, numel(sample_i));
         dhdt = dh./dt_mtx;
         % Normalize the time series
         %thalweg_sample_ht = thalweg_sample_ht./thalweg_sample_ht(1,:);
         subplot(3,9, md_count)
         plot(time(1:end-1), dhdt)
         ylim([-300, 0])
         subplot_title = strrep(folder_dir(i).name, '_',', ');
         title(subplot_title)
%          plot(time, thalweg_sample_ht)
%          title(subplot_title)
%          ylim([0,8000])
    end
end

%% visualize surface and base profiles, and the end of the transient run
folder_dir = dir([pwd '/' foldername]);

figure;
plot_idx = 0;
md_count = 0;
for i = 1:size(folder_dir,1)
    % skip the irrelevant ones
    if ~strcmp(folder_dir(i).name(1), 'm')
        continue
    else
        md_count = md_count + 1;
        try
            md = load([folder_dir(i).folder '/' folder_dir(i).name '/' model_type]).md;
        catch
            continue
        end
        nt = size(md.results.TransientSolution,2);
        
        % geometry parameters
        Lx = max(md.mesh.x);
        Ly = max(md.mesh.y);
        ds = 250; % spacing, 250 meter
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        
        try
            surface = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                           md.results.TransientSolution(end).Surface,...
                                           x, y, NaN);
            base = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                        md.results.TransientSolution(end).Base,...
                                        x, y, NaN);
            bed = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                        md.geometry.bed,...
                                        x, y, NaN);
        catch
            continue
        end
        surface_profile  = surface(mid_i,:);
        base_profile = base(mid_i,:);
        bed_profile = bed(mid_i,:);
        % plot
        subplot(3, 9, md_count)
        plot(x, surface_profile, 'red');hold on;
        plot(x, base_profile,'blue');hold on;
        plot(x, bed_profile, 'black');hold off
        ylim([-600,2000])
        subplot_title = strrep(folder_dir(i).name, '_',', ');
        title(subplot_title)
    end
end