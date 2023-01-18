% In this script, we analyze the force balance of two selected glaciers
% we plot (in GIFs) the force balance of two models side by side
md1_name = "model_W5000_GL0_FC120000";
md2_name = "model_W5000_GL0_FC30000";
model_type = "MISMIP_yangTransient_Calving_MassUnloading.mat";
md1 = load("long_models_yang/" + md1_name + "/" + model_type).md;
md2 = load("long_models_yang/" + md2_name + "/" + model_type).md;
mds = [md1, md2];

gif('plots/force_balance_planview.gif')
figure;
for md_i = 1:length(mds)
    md = mds(md_i);
    index = md.mesh.elements;
    %compute nodal functions coefficients N(x,y)=alpha x + beta y +gamma
    [alpha, beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);
    summation=[1;1;1];
    
    nt = size(md.results.TransientSolution,2);
    Lx = max(md.mesh.x);
    Ly = max(md.mesh.y);
    ds = 50;
    x = 0:ds:Lx;
    y = 0:ds:Ly;
    [X,~] = meshgrid(x, y);
    if rem(size(X,1), 2) == 0
        mid_i = size(X,1)/2;
    else
        mid_i = (size(X,1)+1)/2;
    end
    thalweg_x = X(mid_i,:);
    
    color_length = nt;
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
        linspace(red(2),sth(2),color_length)',...
        linspace(red(3),sth(3),color_length)'];
    fb_ratio_last = 0;
    iter_count = 0;
    % iterate over time
    for i = 30:5:nt
        iter_count = iter_count + 1;
        tauxx = md.results.TransientSolution(i).DeviatoricStressxx;
        tauxy = md.results.TransientSolution(i).DeviatoricStressxy;
        tauxxlist=tauxx(index);
        tauxylist=tauxy(index);
        % get H from vertices to elements
        H = md.results.TransientSolution(i).Thickness;
        H_list = H(index);
        H_list = mean(H_list,2);
        % find directional derivative along x, y
        dtauxxdx=(tauxxlist.*H_list.*alpha)*summation;
        dtauxxdy=(tauxxlist.*H_list.*beta)*summation;
        dtauxydx=(tauxylist.*H_list.*alpha)*summation;
        dtauxydy=(tauxylist.*H_list.*beta)*summation;
        % basal stress; get onto elements
        if size(md.friction.C,2) == 1
            % no sliding law coefficient change
            bs = md.friction.C.^2.*md.results.TransientSolution(i).Vel/md.constants.yts;
        else
            % mass unloading experiment
            bs = md.friction.C(1:end-1,i).^2.*md.results.TransientSolution(i).Vel/md.constants.yts;
        end
        bs_list = bs(index);
        bs = mean(bs_list,2);
        % driving stress
        ds = drivingstress_from_results(md, i);
        
        time = md.results.TransientSolution(i).time;
        plot_title = [md.miscellaneous.name, ', time = ', num2str(time)];
        mask = md.results.TransientSolution(i).MaskOceanLevelset;
        mask = mean(mask(index),2);
    
        % force balance: longitudinal-lateral stress gradient / driving
        % stress
        % The part of d(tauxx)/dx that contributes to the driving stress,
        % we single it out, remove from resistive calculation, and add to
        % the driving stress
        ds_dtauxxdx = zeros(size(dtauxxdx));
        ds_dtauxxdx(dtauxxdx>0) = dtauxxdx(dtauxxdx>0);
        dtauxxdx(dtauxxdx>0) = 0;
        resistive = -(dtauxxdx + dtauxydy);
        fb_ratio = resistive./(ds + ds_dtauxxdx);
        
        fb_ratio(mask<0) = nan;
        
%        subplot(1,length(mds),md_i)

        Lx = max(md.mesh.x);
        Ly = max(md.mesh.y);
        ds = 250;
        x = 0:ds:Lx;
        y = 0:ds:Ly;
        [X,~] = meshgrid(x, y);
        if rem(size(X,1), 2) == 0
            mid_i = size(X,1)/2;
        else
            mid_i = (size(X,1)+1)/2;
        end
        thalweg_x = X(mid_i,:);
        
        mask = md.results.TransientSolution(i).MaskOceanLevelset;
        field_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                            fb_ratio, x, y, NaN);
        field_profile = field_grid(mid_i,:);
        field_prof_smooth = smooth(field_profile(~isnan(field_profile)),15);
        field_profile(~isnan(field_profile)) = field_prof_smooth;
%         delta_field_profile = field_profile - field_profile_last;
%         % update field_profile_last
%         field_profile_last = field_profile;
%         % look for only 15 km behind the grounding line
        if iter_count == 1
            continue 
        else

%             gl_x = locate_groundingline(md,mask);
%             x_keep = find(thalweg_x > gl_x-1.5e4 & thalweg_x < gl_x);
%             % truncate
%             field_profile = field_profile(x_keep);
%             thalweg_x = thalweg_x(x_keep);
%             thalweg_x = thalweg_x - max(thalweg_x);
%         % plot
%             plot(thalweg_x, field_profile, Color=colors_p(i,:));hold on;
%             title(plot_title)
%             ylim([0,1])
%             pause(0.1)
            gif
            plotmodel(md,'data',fb_ratio - fb_ratio_last,'caxis',[0,0.2],'mask',mask)
            fb_ratio_last = fb_ratio;
            pause(0.2)
        end
        %plotmodel(md,'data',-dtauxxdx,'caxis',[0,2e4])
    end
end

