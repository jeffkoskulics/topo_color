function plot_surface(surface)
% Plot surface
surface(:,:,3) = surface(:,:,3) - mean(reshape(surface(:,:,3),...
    numel(surface(:,:,3)),1));

figure;
surf(surface(:,:,1),-surface(:,:,2),surface(:,:,3));
daspect([1 1 1]);
shading interp
colormap([0 0 0]);
light('Position',[-0.2 1 0.6],'Style','Infinite','Color',0.75*[1 1 0.8]);
set(gca,'CameraPosition',[0.25 -1 0.4]);
set(gca,'Projection','perspective');
set(gca,'CameraViewAngle',16);
width = 1920;
height = width / 1.7778;
set(gcf,'OuterPosition',[-10 -10 width height]);
oversize = 1.5;
set(gca,'OuterPosition',[-0.23 -0.2 1.5 1.5]);
set(gca,'XLim',[1.001*surface(1,1,1) surface(1,1280,1)]);
set(gca,'YLim',[surface(1,1,2) surface(1024,1,2)]);
set(gca,'ZLim',[-0.01 0.01]);
set(gcf, 'PaperPositionMode','auto');
set(gca,'FontSize',20);
