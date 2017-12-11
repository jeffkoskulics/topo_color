%Script to plot renderings of the surfaces
directory = 'F:\SIT_December_2011_Results\Wind16\';
list = dir([directory '*.mat']);
N = length(list);

for i = 1:N
    surface = xy_position;
    load([directory list(i).name]);
    surface(:,:,3) = double(elevation);
    plot_surface(double(surface));
    name = [directory 'images\' list(i).name '.tif'];
    print(gcf,'-dtiff',name);
    close(gcf)
end