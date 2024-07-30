function plot_output2D(md, var, varargin)
% Visualize time-dependent variables in 2D 
% Examples: Vel, Elevation, tau_b, dhdt, dudt
%   e.g., plot_output2D(md, 'dhdt')
    nt = size(md.results.TransientSolution,2);
    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)

    ds = 50;
    smooth_window = 10;
    pause_dt = 0.2;

    figure('Position',[100,100, 500,500]);
    for i = 1:nt
            % velocity field
        year = md.results.TransientSolution(i).time;
        plot_title = [md.miscellaneous.name, ', year = ', num2str(year)];
        if strcmp('Vel', var)
            if length(varargin) == 1
                plotmodel(md,'data', md.results.TransientSolution(i).Vel,'title',plot_title,'caxis', varargin{1})
            else
                plotmodel(md,'data', md.results.TransientSolution(i).Vel,'title',plot_title)
            end

        % thickness
        elseif strcmp('Elevation', var)
            if length(varargin) == 1
                plotmodel(md,'data', md.results.TransientSolution(i).Surface-md.results.TransientSolution(1).Surface,'title',plot_title,'caxis',varargin{1})
            else
                plotmodel(md,'data', md.results.TransientSolution(i).Surface-md.results.TransientSolution(1).Surface,'title',plot_title,'caxis',[-300,0],'mask',md.results.TransientSolution(i).MaskIceLevelset<0)
            end

        % basal drag
        elseif strcmp('tau_b', var)
            if ~length(varargin) == 1
                taub_clim = [0,3e5];
            end
            if size(md.friction.C,2) == size(md.results.TransientSolution,2) % not localized basal perturb
                tau_b = md.results.TransientSolution(i).Vel/md.constants.yts.*md.friction.C(1:end-1,i).^2;
            else
                tau_b = md.results.TransientSolution(i).Vel/md.constants.yts.*md.friction.C(1:end-1,(i-1)*10+1).^2;
            end
            plotmodel(md,'data', tau_b,'title',plot_title,'caxis',taub_clim)

        % thickness change
        elseif strcmp('dhdt', var)
            threshold = 0.2; % threshold for dh/dt steady-state; default 0.01 m/a
            if i == 1 % first time step
                plotmodel(md,'data',zeros(size(md.results.TransientSolution(i).Thickness)))
                colorbar
            else
                dh = md.results.TransientSolution(i).Surface - md.results.TransientSolution(i-1).Surface;
                dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
                dhdt = dh/dt;
                dhdt(abs(dhdt) < threshold) = nan;
                time = md.results.TransientSolution(i-1).time;
                title = ['time = ', num2str(time),', threshold at ', num2str(threshold)];
                
                plotmodel(md,'data',dhdt, 'caxis', [-20,10],'title',title)
                pause(0.3)
                colorbar
                %plotmodel(md,'data',dhdt, 'caxis', [-20,20],'title',title,'mask',md.results.TransientSolution(i-1).MaskOceanLevelset>0)
            end

            % velocity change
        elseif strcmp('dudt', var)
            climits = 200;
            if length(varargin) == 1
                threshold = varargin{1};
            else
                threshold = 0.01; % threshold for steady-state
            end
            if i == 1
                plotmodel(md,'data',zeros(size(md.results.TransientSolution(i).Vel)))
            else
                du = md.results.TransientSolution(i).Vel - md.results.TransientSolution(i-1).Vel;
                dt = md.results.TransientSolution(i).time - md.results.TransientSolution(i-1).time;
                dudt = du/dt;
                dudt(abs(dudt) < threshold) = nan;
                time = md.results.TransientSolution(i-1).time;
                title = ['time = ', num2str(time),', threshold at ', num2str(threshold)];
                
                plotmodel(md,'data',dudt, 'caxis', [-climits, climits],'title',title)
            end

        elseif strcmp('ee', var) % effective strain rate
            if length(varargin) == 1
                plotmodel(md,'data', md.results.TransientSolution(i).StrainRateeffective,'title',md.miscellaneous.name,'caxis', varargin{1})
            else
                plotmodel(md,'data', md.results.TransientSolution(i).StrainRateeffective,'title',md.miscellaneous.name)
            end

        elseif strcmp('tau_xx', var) % deviatoric stress
            if length(varargin) == 1
                plotmodel(md,'data', md.results.TransientSolution(i).DeviatoricStressxx,'title',md.miscellaneous.name,'caxis', varargin{1})
            else
                plotmodel(md,'data', md.results.TransientSolution(i).DeviatoricStressxx,'title',md.miscellaneous.name)
                mesh_data = md.results.TransientSolution(i).DeviatoricStressxx;
                grid_data = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
                imagesc(grid_data)
            end

        elseif strcmp('lsg', var) % longitudinal stress gradient
            if length(varargin) == 1
                mesh_data = md.results.TransientSolution(i).DeviatoricStressxx;
                grid_data_tauxx = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
                % H
                mesh_data = md.results.TransientSolution(i).Thickness;
                grid_data_H = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
                Htauxx = grid_data_tauxx.*grid_data_H;
                Htauxx = smooth_2D(Htauxx, 20);
                [lsg, ~] = gradient(Htauxx, 50);
                imagesc(lsg)
                colorbar
                clim(varargin{1})
            else
                mesh_data = md.results.TransientSolution(i).DeviatoricStressxx;
                grid_data_tauxx = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
                % H
                mesh_data = md.results.TransientSolution(i).Thickness;
                grid_data_H = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
                Htauxx = grid_data_tauxx.*grid_data_H;
                Htauxx = smooth_2D(Htauxx, smooth_window);
                [lsg, ~] = gradient(Htauxx, ds);
                imagesc(lsg)
            end

        elseif strcmp('fb', var)
            % in this force balance measure, we show how much driving
            % stress is taken by the longitudianal stress gradient
            mesh_data = md.results.TransientSolution(i).DeviatoricStressxx;
            grid_data_tauxx = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
            % H
            mesh_data = md.results.TransientSolution(i).Thickness;
            grid_data_H = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
            Htauxx = grid_data_tauxx.*grid_data_H;
            Htauxx = smooth_2D(Htauxx, smooth_window);
            [lsg, ~] = gradient(Htauxx, ds);
            % then calculate the driving stress
            mesh_data = md.results.TransientSolution(i).Surface;
            grid_data_s = mesh_to_grid(md.mesh.elements, md.mesh.x, md.mesh.y, mesh_data);
            grid_data_s = smooth_2D(grid_data_s, smooth_window);
            grid_data_alpha = gradient(grid_data_s, ds);
            tau_d = md.materials.rho_ice*md.constants.g.*grid_data_H.*grid_data_alpha;
            % ratio
            fb_ratio = -lsg./tau_d;
            imagesc(fb_ratio)
            clim([-1,1])
            colorbar           
        end
        pause(pause_dt)
        
    end
    hold off
end

