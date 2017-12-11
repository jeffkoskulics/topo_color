%[a,b,row_curv,col_curv] = connection(1024,1280);
system_constants
rows = 1025;
cols = 1281;

%Fourier Plot
    subplot(3,4,[1 2 3 5 6 7 9 10 11]);
    set(gcf,'renderer','ope')
    set(gcf,'OuterPosition',[0 0 1000 1000]);
    %Tukey windowing function
    col_window = tukeywin(cols,0.2);
    row_window = tukeywin(rows,0.2*cols/rows);
    window(:,:,1) = repmat(col_window',rows,1);
    window(:,:,2) = repmat(row_window,1,cols);
    window = window(:,:,1).*window(:,:,2);

    %Window function applied to data
    mean_elevation = mean(mean(elevation));
    g = double(window.*(elevation-mean_elevation));
    L = fft2(g);
    L_shift = fftshift(L);

    %Wavenumber matrix
    k_x = (0:(cols-1)/2);
    k_x((cols-1)/2+2:cols) = -1*((cols-1)/2:-1:1);
    k_x = 2*pi*k_x/(cols*fp_res);
    k_y = (0:(rows-1)/2);
    k_y((rows-1)/2+2:rows) = -1*((rows-1)/2:-1:1);
    k_y = 2*pi*k_y/(rows*fp_res);
    k_y = k_y';
    k_x = repmat(k_x,rows,1);
    k_x = fftshift(k_x);
    k_y = repmat(k_y,1,cols);
    k_y = fftshift(k_y);
    K = (k_x.^2 + k_y.^2).^(1/2);
    slopes = K.^2 .* L_shift;

    Hue = (angle(L_shift) + pi)/(2*pi);
    Value_amplitude = abs(L_shift) / max(reshape(abs(L_shift),numel(L_shift),1));
    Value_slopes = abs(slopes) / max(reshape(abs(slopes),numel(slopes),1));
    H = Hue;
    Value_slopes = 2*Value_slopes;
    H(:,:,2) = Value_slopes;
    H(:,:,2) = H(:,:,2) - (Value_slopes - 1).*double(gt(H(:,:,2),1));
    H(:,:,3) = 1;
    RGB = hsv2rgb(H);
    
    imagesc(RGB)
    title('Hue = Phase, Saturation = Slope amplitude','FontWeight','bold');

    %Tick markers
    x_tick_num = 45;
    y_tick_num = 45;
    x_tick_locations = linspace(1,1280,x_tick_num);
    x_tick_labels = round(k_x(1,uint16(x_tick_locations)));
    y_tick_locations = linspace(1,1024,y_tick_num);
    y_tick_labels = round(k_y(uint16(y_tick_locations),1));
    set(gca,'XTick',x_tick_locations+0.5);
    set(gca,'XTickLabel',x_tick_labels,'FontWeight','bold');
    set(gca,'YTick',y_tick_locations+0.5);
    set(gca,'YTickLabel',y_tick_labels,'FontWeight','bold');
    set(gca,'XLim',[0 160]+640);
    set(gca,'YLim',[-80 80]+512);
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'FontWeight','bold');
    
    %Surface Rendering
    subplot(3,4,8);
    surf(surface(:,:,1),-surface(:,:,2),surface(:,:,3));
    daspect([1 1 1]);
    shading interp
    colormap([0 0 0]);
    light('Position',[-0.2 1 0.6],'Style','Infinite','Color',0.75*[1 1 0.8]);
    set(gca,'CameraPosition',[0.25 -1 0.4]);
    set(gca,'Projection','perspective');
    set(gca,'CameraViewAngle',16);
    set(gca,'XLim',[1.001*surface(1,1,1) surface(1,1280,1)]);
    set(gca,'YLim',[surface(1,1,2) surface(1024,1,2)]);
    set(gca,'ZLim',[-0.01 0.01]);
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    set(gca,'ZTick',[])
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.44 -.30 .2 .2]);
    set(gca,'color', 'none');
    set(gca,'xcolor','white');
    set(gca,'ycolor','white');
    set(gca,'zcolor','white');
    set(gcf, 'PaperPositionMode', 'auto')
    
    