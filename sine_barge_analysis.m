%Sine barge analysis
%Identify uncorrupted pixels 
top    = 99; %top row
bottom = 557; %bottom row
left   = 105; %left column
right  = 966; %right column

%Find indices of uncorrupted and corrupted pixels
good_rows = top:bottom;
good_cols = left:right;
good_pixels = zeros(1025,1281);
good_pixels(good_rows,good_cols) = 1;
good_pixels = gt(good_pixels,0.9);
bad_pixels = ~good_pixels;

ideal = 0.005*sin(2*pi/0.05*xy_position(:,:,1)+0.65);
ideal(reshape(bad_pixels,numel(bad_pixels),1)) = 0;
error = elevation - ideal;

blue = [100*ones(1,100) 100:-1:1];
red  = [1:1:100 100*ones(1,100)];
green = [1:1:100 100:-1:1];

RGB = ([red' green' blue']);

figure; imagesc(error,[-0.005 0.005])
colormap(RGB/100)

figure; plot(xy_position(325,:,1),elevation(325,:)); 
hold on; plot(xy_position(325,:,1),ideal(325,:),'red')