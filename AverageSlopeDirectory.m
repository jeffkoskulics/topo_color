%Slope histograms for a directory

directory = 'C:\Users\Jeff\Documents\A Research\WallopsDec2011\Wind16\Retrieval\Wind16\';
list = dir([directory '*.mat']);
load([directory list(length(list)).name]);
z = xy_position;
 %Allocate space for the histogram
 S = zeros(length(y_bins),length(x_bins));

%function [s sx sy] = slope_histogram(z)
for k = 1:1000
    load([directory list(k).name]);
    z(:,:,3) = elevation;
    %compute surface normals
    [u_sn, r_sn] = surface_normals(z);

    %x-component of slope
     sx = reshape(u_sn(:,:,1)./u_sn(:,:,3),numel(u_sn(:,:,1)),1);
    %y-component of slope
     sy = reshape(u_sn(:,:,2)./u_sn(:,:,3),numel(u_sn(:,:,1)),1);
    %slope magnitude
     s = (sx.^2 +sy.^2).^0.5;

     %Set up histogram bins
     max_x_slope = 0.5774;
     min_x_slope = -0.5774;
     num_x_bins = 100;
     max_y_slope = 0.5774;
     min_y_slope = -0.5774;
     num_y_bins = 100;
     x_bins = linspace(min_x_slope,max_x_slope,num_x_bins);
     y_bins = linspace(min_y_slope,max_y_slope,num_y_bins);

     %Calculate x- and y-bins for all data points
     row_index = round((sy-y_bins(1))/(y_bins(2)-y_bins(1)))+1;
     col_index = round((sx-x_bins(1))/(x_bins(2)-x_bins(1)))+1;

     %Add each data point to its appropriate histogram bin
     for i = 1:length(row_index)
         if ((0 < row_index(i) && row_index(i) <= num_y_bins) && ...
                 (0 < col_index(i) && col_index(i) <= num_x_bins))
            S(row_index(i),col_index(i)) = S(row_index(i),col_index(i)) +1; 
         end
     end
     disp(k)
end