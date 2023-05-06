function [fig, ax] = plot_receiver(receiver, sim_reg, S, w0_vect, r_0_L)
fig = figure();
for conta=1:length(S)
    S1=S(conta);
    [U_ric] = receiver.U;
    
    ax = subplot(1,1,conta);

    hold on;
    imagesc(sim_reg.x*1e2, sim_reg.y*1e2, U_ric);
    viscircles([0 0],r_0_L/2*1e2,'Color','k','linestyle','--','linewidth',1, 'Color', 'r');
    text(r_0_L/2, r_0_L/2*1e2, 'r0', 'Color', 'w');
    hold off;

    ratio=2*w0_vect(1)/r_0_L;
    title(['D/r_0 = ', num2str(ratio)]);
    % set(gca, 'Ydir', 'normal');
    xlabel('x [cm]');
    ylabel('y [cm]');
    set(gca,'fontsize',14,'linewidth',1);
    axis("tight")
    axis("square")
    box("on")
end
end