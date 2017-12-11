%fft analysis
[rows cols] = size(elevation);
fp_res = 0.000297;

%Tukey windowing function
col_window = tukeywin(cols,0.2);
row_window = tukeywin(rows,0.2*cols/rows);
window(:,:,1) = repmat(col_window',rows,1);
window(:,:,2) = repmat(row_window,1,cols);
window = window(:,:,1).*window(:,:,2);

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

%{ 
Plot the windowed surface
surface = xy_position;
surface(:,:,3) = window(:,:,3).*elevation;
plot_surface(surface);
%}

%Window function applied to data
mean_elevation = mean(mean(elevation));
z = window.*(elevation-mean_elevation);
L = fft2(z);
L_shift = fftshift(L);
Hue = (angle(L_shift) + pi)/(2*pi);
Value_amplitude = abs(L_shift) / max(reshape(abs(L_shift),numel(L_shift),1));
Value_slopes = abs(slopes) / max(reshape(abs(slopes),numel(slopes),1));
H = Hue;
Value_slopes = 2*Value_slopes;
H(:,:,2) = Value_slopes;
H(:,:,2) = H(:,:,2) - (Value_slopes - 1).*double(gt(H(:,:,2),1));
H(:,:,3) = 1;
RGB = hsv2rgb(H);
figure; imagesc(RGB);
title('Hue = Phase, Saturation = Slope amplitude');

%Tick markers
x_tick_num = 45;
y_tick_num = 45;
x_tick_locations = linspace(1,1280,x_tick_num);
x_tick_labels = round(k_x(1,uint16(x_tick_locations)));
y_tick_locations = linspace(1,1024,y_tick_num);
y_tick_labels = round(k_y(uint16(y_tick_locations),1));
set(gca,'XTick',x_tick_locations+0.5);
set(gca,'XTickLabel',x_tick_labels);
set(gca,'YTick',y_tick_locations+0.5);
set(gca,'YTickLabel',y_tick_labels);
set(gca,'XLim',[0 160]+640);
set(gca,'YLim',[-80 80]+512);
set(gca,'DataAspectRatio',[1 1 1]);

%{
slopes = K.^2 .* L_shift;

H(:,:,3) = Value_slopes;
RGB1 = hsv2rgb(H);
figure; imagesc(RGB1);
title('Hue = Phase, Value = Slope amplitude');

%}