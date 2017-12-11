function plot_surface_diff(surface,previous_surface)

%Calc the diff
surface = cat(3,surface(:,:,1:2),surface(:,:,3)-previous_surface(:,:,3));

% Plot surface
%surface = reshape(surface,size(surface,1)*size(surface,2),3);
bad_pixels = lt(surface(:,:,3),-.5) | gt(surface(:,:,3),.5);
temp = surface(:,:,3);
surface_fixed = zeros(size(surface,1),size(surface,2),1);
surface_fixed(bad_pixels) = NaN;
surface_fixed(~bad_pixels) = temp(~bad_pixels);

figure;
surf(surface(:,:,1),surface(:,:,2),surface_fixed(:,:));
daspect([1 1 1]);
shading interp
colormap([0 0 0]);
light('Position',[-0.2 1 0.6],'Style','Infinite','Color',0.75*[1 1 0.8]);
set(gca,'CameraPosition',[0 -1 0.45]);
set(gca,'Projection','perspective');
set(gca,'CameraViewAngle',18);
width = 1920;
height = width / 1.7778;
set(gcf,'OuterPosition',[-10 -10 width height]);
oversize = 1.5;
set(gca,'OuterPosition',[-0.27 -0.1 1.5 1.5]);
set(gca,'ZLim',0+[-.05 .05]);
set(gcf, 'PaperPositionMode','auto');
%set(gcf,'Position',get(gcf,'OuterPosition')*oversize); 
%set(gcf,'Renderer','zbuffer')
%set(gca,'FaceLighting','phong','AmbientStrength',0.5);
%set(h,'AmbientLightColor',[0.2 0.2 0.8]);