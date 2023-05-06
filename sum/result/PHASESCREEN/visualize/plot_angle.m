function [fig, ax, cbar] = plot_angle(U, sim_reg, NameValueargs)
    arguments
        U (:, :) {mustBeNumeric};
        sim_reg SimulationRegion;
        NameValueargs.unwrap {mustBeNumericOrLogical} = false;
    end

    phase = angle(U);
    % Normalize to the peak intensity
    if(NameValueargs.unwrap)

    else

    fig = figure();
    ax = subplot(1, 1, 1);
    imagesc(ax, sim_reg.x*1e3, sim_reg.y*1e3, phase);
    title('Phase front at Receiver');
    xlabel('x [mm]');
    ylabel('y [mm]');
    axis(ax, "square");
    colormap(ax, "hot");
    cbar = colorbar(ax);
    cbar.Limits = [-pi, pi];
    cbar.Ticks = -pi:pi/2:pi;
    cbar.TickLabels = {'-\pi', '-\pi/2','0','\pi/2','\pi'};
end

