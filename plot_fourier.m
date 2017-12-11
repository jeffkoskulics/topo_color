function plot_fourier(elevation,xy_position)

z = xy_position;
z(:,:,3) = elevation;
elevation = z;

%Fit a least squares plane
x = reshape(elevation(:,:,1),numel(elevation(:,:,1)),1);
y = reshape(elevation(:,:,2),numel(elevation(:,:,1)),1);
lin_sys = [x y ones(length(x),1)];
a = lin_sys\reshape(elevation(:,:,3),numel(elevation(:,:,1)),1);

%Subtract least squares plane from data
elevation(:,:,3) = elevation(:,:,3) - a(1)*elevation(:,:,1) - ...
    a(2)*elevation(:,:,2) - a(3);

%plot_surface(elevation);
%plot_surface(cat(3,elevation(:,:,1:2),a(1)*elevation(:,:,1) + ...
%    a(2)*elevation(:,:,2) + a(3)));
elevation = elevation(:,:,3);

%fft analysis
[rows cols] = size(elevation);
fp_res = 0.000297;

%Tukey windowing function
window_width = 0.25;
col_window = tukeywin(cols,window_width);
row_window = tukeywin(rows,window_width*cols/rows);
window(:,:,1) = repmat(col_window',rows,1);
window(:,:,2) = repmat(row_window,1,cols);
window = window(:,:,1).*window(:,:,2);

%Window function applied to data
z = window.*elevation;
z = z - mean(mean(z));
L = fft2(z)/((rows/2)*(cols/2));
FT = fftshift(L);

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
amplitude = abs(FT);
slopes = (atan(K.^2 .* abs(FT))*180/pi); %Should it be abs(FT).^2 ?
figure; plot(K(513,:),slopes(511:514,:))
title('slopes')
figure; plot(K(513,:),amplitude(511:514,:))
title('elevation')

Hue = (angle(FT) + pi)/(2*pi);
Saturation = slopes/(max(max(slopes)))*1;
Saturation = Saturation - (Saturation - 1).*double(gt(Saturation,1));
H = Hue;
H(:,:,2) = Saturation;
H(:,:,3) = 1;
RGB = hsv2rgb(H);

figure;
imagesc(RGB)
title('Hue = Phase, Saturation = Slope amplitude','FontSize',16);

%Tick markers
x_tick_num = 45;
y_tick_num = 45;
%x_tick_locations = linspace(1,1280,x_tick_num);
x_tick_locations = [640 678 714 750 787];
%x_tick_labels = round(k_x(1,uint16(x_tick_locations)));
x_tick_labels = [0 600 1200 1800 2400];
%y_tick_locations = linspace(1,1024,y_tick_num);
y_tick_locations = [435 474 513 552 591];
%y_tick_labels = round(k_y(uint16(y_tick_locations),1)/100)*100;
y_tick_labels = [1600 800 0 -800 -1600];
set(gca,'XTick',x_tick_locations+0.5);
set(gca,'XTickLabel',x_tick_labels,'FontSize',16);
xlabel('upwind wavenumber (rad m^-^1)','FontSize',16)
set(gca,'YTick',y_tick_locations+0.5);
set(gca,'YTickLabel',y_tick_labels);
ylabel('crosswind wavenumber (rad m^-^1)','FontSize',16)
%set(gca,'XLim',[0 160]+640);

%set(gca,'YLim',[-80 80]+512);
set(gca,'DataAspectRatio',[1 1 1]);



%{
slopes = K.^2 .* FT;
H(:,:,3) = Value_slopes;
RGB1 = hsv2rgb(H);
figure; imagesc(RGB1);
title('Hue = Phase, Value = Slope amplitude');

%}