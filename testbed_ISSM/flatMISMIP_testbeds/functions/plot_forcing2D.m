function plot_forcing2D(md, var, varargin)
%PLOT_FORCING2D make animated plots of forcing fields
    warning_id = 'MATLAB:handle_graphics:exceptions:SceneNode';
    warning('off',warning_id)
    if strcmp('B', var)
        nt = size(md.materials.rheology_B, 2);
        figure;
        for i = 1:nt
            % velocity field
            year = md.materials.rheology_B(end,i);
            plot_title = [md.miscellaneous.name, ', year = ', num2str(year)];
            if length(varargin) == 1
                plotmodel(md,'data', md.materials.rheology_B(1:end-1,i),'title',plot_title, 'caxis', varargin{1})
            else
                plotmodel(md,'data', md.materials.rheology_B(1:end-1,i),'title',plot_title)
            end
        end
    elseif strcmp('B_ratio', var)
        nt = size(md.materials.rheology_B, 2);
        figure;
        B0 = md.materials.rheology_B(1:end-1,1);
        for i = 1:nt
            % velocity field
            year = md.materials.rheology_B(end,i);
            plot_title = [md.miscellaneous.name, ', year = ', num2str(year)];
            plotmodel(md,'data', md.materials.rheology_B(1:end-1,i)./B0,'title',plot_title,'caxis',[0,1])
        end

    elseif strcmp('C_ratio', var)
        nt = size(md.friction.C, 2);
        figure;
        C0 = real(md.friction.C(1:end-1,1));
        for i = 1:nt
            % velocity field
            year = md.friction.C(end,i);
            plot_title = [md.miscellaneous.name, ', year = ', num2str(year)];
            plotmodel(md,'data', real(md.friction.C(1:end-1,i))./C0,'title',plot_title,'caxis',[0,1])
        end

    end
end

