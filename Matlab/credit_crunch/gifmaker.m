
load huggett_transition_creditcrunch.mat
for n=1:length(gg)/8 % First 50 periods
    plot(a,gg{n}(1:I),'LineWidth',3,'color',[0 0 153/255])
    hold on
    plot(a,gg{n}(I+1:2*I),'LineWidth',3,'color',[1 0 0])
    hold off
    ylim([0 3.5])
    xlim([amin-num*da 0.5])
    title(sprintf('test (t=%0.1f)',n));
    xlabel('$Wealth$','interpreter','latex')
    legend('g_1(a,t)','g_2(a,t)','interpreter','latex')
    legend('boxoff')
    alpha(0.15)
    grid on
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    outfile = 'creditcrunch.gif';
        if n==1
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
    end
end


