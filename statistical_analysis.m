%Statistical analysis
[a b c d] = connection(1024,1280);
fp_res = xy_position(1,2,1) - xy_position(1,1,1);

%Finite difference matrices for slopes and curvature
x_slopes = (a * fp_res^-1)*180/pi;
y_slopes = (b * fp_res^-1)*180/pi;
x_curves = c * fp_res^-2;
y_curves = d * fp_res^-2;

%Topography vector
z = double(reshape(elevation,1025*1281,1));

x_s = x_slopes * z;
y_s = y_slopes * z;
x_c = x_curves * z;
y_c = y_curves * z;

elev_bins  = linspace(min(z),max(z),200);
elev_bin_width = elev_bins(2) - elev_bins(1);
slope_bins = linspace(min([x_s;y_s]),max([x_s;y_s]),100);
slope_bin_width = slope_bins(2) - slope_bins(1);
curve_bins = linspace(min([x_c;y_c]),max([x_c;y_c]),300);
curve_bin_width = curve_bins(2) - curve_bins(1);

elev_hist = hist(z,elev_bins);
elev_hist = elev_hist / (sum(elev_hist)*elev_bin_width);
figure; plot(elev_bins,elev_hist);

x_slope_hist = hist(x_s,slope_bins);
x_slope_hist = x_slope_hist / (sum(x_slope_hist)*slope_bin_width);
y_slope_hist = hist(y_s,slope_bins);
y_slope_hist = y_slope_hist / (sum(y_slope_hist)*slope_bin_width);
figure; plot(slope_bins,x_slope_hist,'red')
hold on
plot(slope_bins,y_slope_hist,'green')

x_curve_hist = hist(x_c,curve_bins);
x_curve_hist = x_curve_hist / (sum(x_curve_hist)*curve_bin_width);
y_curve_hist = hist(y_c,curve_bins);
y_curve_hist = y_curve_hist / (sum(y_curve_hist)*curve_bin_width);
figure; plot(curve_bins,x_curve_hist, 'red');
hold on
plot(curve_bins,y_curve_hist,'green');

disp(['elev mean ' num2str(mean(z)) ' std ' num2str(std(z))])
disp(['x slopes ' num2str(mean(x_s)) ' std ' num2str(std(x_s))])
disp(['y slopes ' num2str(mean(y_s)) ' std ' num2str(std(y_s))])
disp(['x curves ' num2str(mean(x_c)) ' std ' num2str(std(x_c))])
disp(['y curves ' num2str(mean(y_c)) ' std ' num2str(std(y_c))])

surface = xy_position;
surface(:,:,3) = elevation;
plot_surface(surface)
figure; plot(xy_position(512,:,1),elevation(512,:))