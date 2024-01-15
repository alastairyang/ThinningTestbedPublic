function  export_graph(md, model_index, model_type)
%EXPORT_GRAPH It exports and saves requested plots so that I can monitor
%the simulation results as it goes
%
%   Input:
%       md: ISSM model structure
%       model_index [str]: ...
%       model_type  [str]: "ss","spinup","t"

    figsave_path = '/Users/donglaiyang/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alastair’s_MacBook_Pro/Buffalo/Research/git_research/testbed_ISSM/my_testbeds/Graphs/ensemble_monitor/';
    if ~isa(model_index, 'char')
        model_index = convertStringsToChars(model_index);
    end
    if ~isa(model_type, 'char')
        model_index = convertStringsToChars(model_type);
    end
    name = ['model_',model_index, '_', model_type];
    nt = size(md.results.TransientSolution,2);
    
    % make plots:
    %       plot 1: the timeseries of mean ice thickness (to monitor if the
    %       relaxation is completed)
    [geometry, ~] = query_data(model_index, model_type);
    syn = testbed_data(geometry{1});
    X = syn.X;
    Y = syn.Y;
    x = syn.x;
    y = syn.y;
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);
    
    % only plot the first and the last timestep
    % we don't want to display the figure -> save it
    f = figure('visible','off','Position',[50,50,1500, 700]);
    interval = max(floor(nt*0.01),1);
    selected_t = 1:interval:nt;
    sample_i = 1:5:100;
    thalweg_sample_ht = [];
    
    subplot(1,2,1)
    for i = selected_t
        thalweg_h_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                            md.results.TransientSolution(i).Surface,...
                                            x, y, NaN);
        thalweg_sample_h  = thalweg_h_grid(mid_i,sample_i);
        thalweg_sample_ht = [thalweg_sample_ht; thalweg_sample_h];

    end
    % Normalize the time series
    thalweg_sample_ht = thalweg_sample_ht./thalweg_sample_ht(1,:);
    plot(selected_t, thalweg_sample_ht)
    title('Surface elevation at various points on the centerline')

    %       plot 2: Ice thickness profile at the first and the last timestep 
    subplot(1,2,2)
    selected_t_startend = [selected_t(1), selected_t(end)];
    labels = ["t = 0", "", "t = final time",""];
    run_count = 0;
    for i = selected_t_startend
        run_count = run_count + 1;
        thalweg_h_grid = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x, md.mesh.y,...
                                            md.results.TransientSolution(i).Surface,...
                                            x, y, NaN);
        thalweg_base_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                            md.results.TransientSolution(i).Surface - md.results.TransientSolution(i).Thickness,...
                                            x, y, NaN);
        thalweg_h_line = thalweg_h_grid(mid_i,:);
        thalweg_base_line = thalweg_base_grid(mid_i,:);
        thalweg_both = [thalweg_h_line; thalweg_base_line];
        plot(thalweg_x, thalweg_both); hold on    
        pause(0.2)
    end
    legend(labels)
    title('Lateral profiles')

    saveas(f, [figsave_path, name,'.png']);

    
end

