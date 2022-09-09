function [] = validationPlot(results,range)
bf = log2(results(:,1));
bs = results(:,2);
error = log2(results(:,13));
vali_fit = fit([bf,bs],error,'Lowess','Normalize','on','Robust','LAR');
plot(vali_fit,'Style','contour')
grid(gca,'off');
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontSize',6);
set(gca,'CLim',range,'Layer','top');
end