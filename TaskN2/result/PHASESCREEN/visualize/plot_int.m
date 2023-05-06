function [fig, ax, cbar] = plot_int(U, sim_reg, NameValueargs)
    arguments
        U (:, :) {mustBeNumeric};
        sim_reg SimulationRegion;
        NameValueargs.normalize {mustBeNumericOrLogical} = true;
    end

    % Normalize to the peak intensity
    if(NameValueargs.normalize)
        img = abs(U).^2/max(abs(U).^2, [], "all");
    else
        img = abs(U).^2;
    end

    fig = figure();
    ax = subplot(1, 1, 1);
    imagesc(ax, sim_reg.x*1e3, sim_reg.y*1e3, img);
    title('Captured intensity at Receiver');
    xlabel('x [mm]');
    ylabel('y [mm]');
    axis(ax, "square");
    colormap(ax, "hot");
    cbar = colorbar(ax);
    cbar.Label.String = "$I$";
    cbar.Label.Interpreter = "latex";
end

