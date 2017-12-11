function RetrievalAnalysis(parent_dir,data_cases)
%This script provides a series of analysis plots on the retrieval surfaces
% and the data images.
%A specific directory structure is required:
% parent_dir
%   - Data_Case_0
%       - image000.tif
%       - image001.tif
%       - ...
%   - Data_Case_1
%       - image000.tif
%       - image001.tif
%       - ...
%   - Retrieval
%       - Data_Case_0
%           - image000.mat
%           - image001.mat
%           - ...
%       - Data_Case_1
%           - image000.mat
%           - image001.mat
%           - ...
%   - flat_model_image.mat (the flat surface model image used during analysis)
%
% Format the input as follows:
% parent_dir: The directory that holds case data directories and the
%             "Retrieval" directory output by ElevationBatch
% data_cases: input as a single column cell vector
%             e.g. {'Wind 16';'Wind500000';'N2';...;...}
%             NOTICE: semicolon and curly braces

% Check to see if directory ends with a '\', and convert windows style to
% unix style
[part,remainder] = strtok(parent_dir,'/\');
if isunix
    parent_dir = ['/',part];
else
    parent_dir = part;
end

while ~isempty(part)
    [part,remainder] = strtok(remainder,'/\');
    parent_dir = [parent_dir,'/',part];
end

%Compile the list of directories to analyze
display('Analyzing the following data cases')
for i = 1:size(data_cases,1)
    directory_list{i} = [parent_dir,data_cases{i},'/'];
    display(directory_list{i})
end

%Check for output directory existance
if ~exist([parent_dir,'RetrievalAnalysis/'],'dir')
    display('Output directory does not exist')
    display('Creating output directories')
    mkdir(parent_dir,'RetrievalAnalysis/');
end

%Used for curvature
[a,b,row_curv,col_curv] = connection(1024,1280);
flat_model_image = load([parent_dir,'flat_model_image.mat']);
flat_model_image = flat_model_image.flat_model_image;
fp_res = 0.0002898; %in meters on the focal plane
gain =[];
offset = [];
data = [];
elevation =[];
rows = 1025;
cols = 1281;

clear File_List
for z = 1:size(directory_list,2)
    
    display(['Processing in directory: ',directory_list{z}])
    
    %Check for output directory existance
    sub_dir = strrep(directory_list{z},parent_dir,'');
    if ~exist([parent_dir,'RetrievalAnalysis/',sub_dir],'dir')
        display('Output directory does not exist')
        display('Creating output directory')
        mkdir([parent_dir,'RetrievalAnalysis/',sub_dir]);
    end
    if ~exist([parent_dir,'RetrievalAnalysis/',sub_dir,'Multiplot-Rendering/'],'dir')
        mkdir([parent_dir,'RetrievalAnalysis/',sub_dir,'Multiplot-Rendering/'])
    end
    if ~exist([parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Diff/'],'dir')
        mkdir([parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Diff/'])
    end
    if ~exist([parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Fourier/'],'dir')
        mkdir([parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Fourier/'])
    end
    if ~exist([parent_dir,'RetrievalAnalysis/',sub_dir,'Data-Diff/'],'dir')
        mkdir([parent_dir,'RetrievalAnalysis/',sub_dir,'Data-Diff/'])
    end
    
    %Load correction matrices
    display('Loading new xy_position.mat...')
    load([parent_dir,'Retrieval/',sub_dir,'xy_position.mat'])
    
    display(['Loading new gain.mat and offset.mat...',10])
    load([parent_dir,sub_dir,'gain.mat'])
    load([parent_dir,sub_dir,'offset.mat'])
    
    File_List = dir([parent_dir,'Retrieval/',sub_dir,'/image*.mat']);
    
    for k = 1:length(File_List)
    
        %Load the .mat elevation retrieval and respective .tif data image
        load([parent_dir,'Retrieval/',sub_dir,File_List(k).name])
        surface = cat(3,xy_position,double(elevation));
        
        data = imread([parent_dir,sub_dir,strrep(File_List(k).name,'.mat','.tif')],'tif');
        data = double(demosaic(uint16((double(data) - offset)),'gbrg')).*gain;
        diff_data = data-flat_model_image;
        diff_data = ((diff_data*128)/6000)+128; %scaled into the 8bit range
        
        %Plot of rendering, elevation profile, elevation hist, slope hist, curv hist
        plot_rendering_profile_hists()
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Multiplot-Rendering/',strrep(File_List(k).name,'.mat','-Multiplot-Rendering.eps')],'eps')
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Multiplot-Rendering/',strrep(File_List(k).name,'.mat','-Multiplot-Rendering.jpg')],'jpg')
        close(gcf)
        
        %Plot of rendering and data image
        plot_rendering_data_image()
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Diff/',strrep(File_List(k).name,'.mat','-Rendering-Diff.eps')],'eps')
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Diff/',strrep(File_List(k).name,'.mat','-Rendering-Diff.jpg')],'jpg')
        close(gcf)
        
        %Plot data image and diff image
        plot_data_and_diff()
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Data-Diff/',strrep(File_List(k).name,'.mat','-Data-Diff.eps')],'eps')
        saveas(gcf,[parent_dir,'RetrievalAnalysis/',sub_dir,'Data-Diff/',strrep(File_List(k).name,'.mat','-Data-Diff.jpg')],'jpg')
        close(gcf)
        
        %Plot rendering with fourier
        plot_rendering_fourier()
        print(gcf,'-painters','-dtiff',[parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Fourier/',strrep(File_List(k).name,'.mat','-Rendering-Fourier.tif')])
        %print(gcf,'-painters','-depsc',[parent_dir,'RetrievalAnalysis/',sub_dir,'Rendering-Fourier/',strrep(File_List(k).name,'.mat','-Rendering-Fourier.eps')])
        close(gcf)
        
    end

    display('Processing complete... ending')
end

function plot_rendering_profile_hists()
    %Surface Rendering
    subplot(3,3,[1 2 4 5])
    set(gcf,'OuterPosition',[0 0 1000 600]);
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
    set(gca,'OuterPosition',[.16 .31 .4 1]);
    set(get(gca,'xlabel'),'Units','normalized','Position',[.5 -.03],'String','Upwind (m)','FontWeight','bold')
    set(get(gca,'ylabel'),'Units','normalized','Position',[.9 0],'String','Crosswind (m)','FontWeight','bold')
    set(get(gca,'zlabel'),'Units','normalized','Position',[-.09 .4],'String','Elevation (m)','FontWeight','bold')
    set(get(gca,'title'),'Units','normalized','Position',[-.02 .80],'String',[strrep(data_cases{z},'_','\_'),10,strrep(strrep(File_List(k).name,'.mat',''),'_','\_')],'FontWeight','bold')
    set(gca,'FontWeight','bold');

    %Elevation cross-section through 0m crosswind
    subplot(3,3,[7 8])
    plot(xy_position(512,:,1),elevation(512,:)*1000,'color','black','LineWidth',1.5)
    set(gca,'OuterPosition',[.04 .06 .60 .45]);
    set(gca,'XLim',[-0.2 0.2]);
    set(gca,'YLim',[-10 10]);
    xlabel('Upwind (m)','FontWeight','bold')
    ylabel('Elevation (mm)','FontWeight','bold')
    set(gca,'FontWeight','bold');

    %Elevation Histogram
    subplot(3,3,3)
    bin_width = .1;
    x_limit = [-11 11];
    h = hist(reshape(elevation,size(elevation,1)*size(elevation,2),1)*1000,x_limit(1):bin_width:x_limit(2));
    stairs(x_limit(1):bin_width:x_limit(2),h/(bin_width*numel(elevation)),'color','black','LineWidth',1.5)
    set(gca,'XLim',[x_limit(1)+1 x_limit(2)-1]);
    set(gca,'YLim',[0 .4]);
    elev_mean = mean(mean(elevation));
    elev_std = std(reshape(elevation,size(elevation,1)*size(elevation,2),1));
    xlabel('Elevation (mm)','FontWeight','bold')
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.01 -.01 .08 -.003]);
    set(get(gca,'ylabel'),'Units','normalized','Position',[.19 .5],'String',['Relative',10,'Frequency',10,'(mm^-^1)'],'FontWeight','bold')
    set(get(gca,'title'),'Units','normalized','Position',[.99 .70],'String',...
        ['\mu =',strrep(strrep(num2str(elev_mean,'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\sigma =',strrep(strrep(num2str(elev_std,'%+.2E'),'E+0','E+'),'E-0','E-')]...
        ,'FontWeight','bold','horizontalAlignment','right','FontName','FixedWidth')
    set(gca,'FontWeight','bold');

    %Slope Histogram
    subplot(3,3,6)
    bin_width = .5;
    x_limit = [-26 26];
    u_sn = surface_normals(cat(3,xy_position,double(elevation)));
    slope(:,:,1) = (180/pi)*atan(u_sn(:,:,1)./u_sn(:,:,3));
    slope(:,:,2) = (180/pi)*atan(u_sn(:,:,2)./u_sn(:,:,3));
    s_hist_x = hist(reshape(slope(:,:,1),size(slope,1)*size(slope,2),1),x_limit(1):bin_width:x_limit(2));
    s_hist_y = hist(reshape(slope(:,:,2),size(slope,1)*size(slope,2),1),x_limit(1):bin_width:x_limit(2));
    stairs(x_limit(1):bin_width:x_limit(2),s_hist_x/(bin_width*numel(slope(:,:,1))),'color','black','LineWidth',1.5)
    hold on
    stairs(x_limit(1):bin_width:x_limit(2),s_hist_y/(bin_width*numel(slope(:,:,2))),':k','LineWidth',2.5)
    hold off
    slope_mean(1) = mean(mean(slope(:,:,1)));
    slope_mean(2) = mean(mean(slope(:,:,2)));
    slope_stdev(1) = std(reshape(slope(:,:,1),size(slope,1)*size(slope,2),1));
    slope_stdev(2) = std(reshape(slope(:,:,2),size(slope,1)*size(slope,2),1));
    set(gca,'XLim',[x_limit(1)+1 x_limit(2)-1]);
    set(gca,'YLim',[0 .35]);
    xlabel('Slope (deg)','FontWeight','bold')
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.01 -.01 .08 -.003]);
    set(get(gca,'ylabel'),'Units','normalized','Position',[.19 .5],'String',['Relative',10,'Frequency',10,'(deg^-^1)'],'FontWeight','bold')
    set(get(gca,'title'),'Units','normalized','Position',[.99 .42],'String',...
        ['\mu x=',strrep(strrep(num2str(slope_mean(1),'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\mu y=',strrep(strrep(num2str(slope_mean(2),'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\sigma x=',strrep(strrep(num2str(slope_stdev(1),'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\sigma y=',strrep(strrep(num2str(slope_stdev(2),'%+.2E'),'E+0','E+'),'E-0','E-')]...
        ,'FontWeight','bold','horizontalAlignment','right','FontName','FixedWidth')
    set(gca,'FontWeight','bold');

    %Curvature Histogram
    subplot(3,3,9)
    bin_width = 1;
    x_limit = [-101 101];
    curv_down = (1/fp_res^2)*row_curv*reshape(double(elevation),size(elevation,1)*size(elevation,2),1);
    curv_cross = (1/fp_res^2)*col_curv*reshape(double(elevation),size(elevation,1)*size(elevation,2),1);
    hist_down = hist(reshape(curv_down,size(curv_down,1)*size(curv_down,2),1),x_limit(1):bin_width:x_limit(2));
    hist_cross = hist(reshape(curv_cross,size(curv_cross,1)*size(curv_cross,2),1),x_limit(1):bin_width:x_limit(2));
    stairs(x_limit(1):bin_width:x_limit(2),hist_down/(bin_width*length(curv_down)),'color','black','LineWidth',1.5)
    hold on
    stairs(x_limit(1):bin_width:x_limit(2),hist_cross/(bin_width*length(curv_cross)),':k','LineWidth',2.5)
    hold off
    curv_d_mean = mean(curv_down);
    curv_c_mean = mean(curv_cross);
    curv_d_std = std(curv_down);
    curv_c_std = std(curv_cross);
    set(gca,'XLim',[x_limit(1)+1 x_limit(2)-1]);
    set(gca,'YLim',[0 .08]);
    xlabel('Curvature (m^-^1)','FontWeight','bold')
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.01 -.01 .08 -.003]);
    set(get(gca,'ylabel'),'Units','normalized','Position',[.19 .5],'String',['Relative',10,'Frequency',10,'(m)'],'FontWeight','bold')
    set(get(gca,'title'),'Units','normalized','Position',[.99 .42],'String',...
        ['\mu x=',strrep(strrep(num2str(curv_d_mean,'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\mu y=',strrep(strrep(num2str(curv_c_mean,'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\sigma x=',strrep(strrep(num2str(curv_d_std,'%+.2E'),'E+0','E+'),'E-0','E-'),10,...
        '\sigma y=',strrep(strrep(num2str(curv_c_std,'%+.2E'),'E+0','E+'),'E-0','E-')]...
        ,'FontWeight','bold','horizontalAlignment','right','FontName','FixedWidth')
    set(gca,'FontWeight','bold');

    set(gcf, 'PaperPositionMode', 'auto')

end

function plot_rendering_data_image()
    %Surface Rendering
    subplot(1,3,[1 2]);
    set(gcf,'OuterPosition',[0 0 1000 600]);
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
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.01 0 0 -.1]);
    set(get(gca,'xlabel'),'Units','normalized','Position',[.5 -.03],'String','Upwind (m)','FontWeight','bold')
    set(get(gca,'ylabel'),'Units','normalized','Position',[.9 0],'String','Crosswind (m)','FontWeight','bold')
    set(get(gca,'zlabel'),'Units','normalized','Position',[-.09 .4],'String','Elevation (m)','FontWeight','bold')
    set(get(gca,'title'),'Units','normalized','Position',[-.02 .80],'String',[strrep(data_cases{z},'_','\_'),10,strrep(strrep(File_List(k).name,'.mat',''),'_','\_')],'FontWeight','bold')
    set(gca,'FontWeight','bold');

    %Data Diff Image
    subplot(1,3,3);
    imagesc(uint8(diff_data))
    daspect([1 1 1]);
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.01 -.10 .10 0]);
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    set(gcf, 'PaperPositionMode', 'auto')
end

function plot_data_and_diff()
    %Data Image
    subplot(1,2,1);
    set(gcf,'OuterPosition',[0 0 720 500]);
    set(get(gca,'title'),'Units','normalized','Position',[-.02 .80],'String',[strrep(data_cases{z},'_','\_'),10,strrep(strrep(File_List(k).name,'.mat',''),'_','\_')],'FontWeight','bold')
    imagesc(uint16(data))
    daspect([1 1 1]);
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.13 0 .15 0]);
    set(gca,'YTick',[])
    set(gca,'XTick',[])

    %Data Diff Image
    subplot(1,2,2);
    imagesc(uint8(diff_data))
    daspect([1 1 1]);
    set(gca,'OuterPosition',get(gca,'OuterPosition') + [-.06 0 .15 0]);
    set(gca,'YTick',[])
    set(gca,'XTick',[])
end

function plot_rendering_fourier()
    %Fourier Plot
    subplot(3,4,[1 2 3 5 6 7 9 10 11]);
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
    set(get(gca,'title'),'String',[strrep(data_cases{z},'_','\_'),' | ',strrep(strrep(File_List(k).name,'.mat',''),'_','\_')],'FontWeight','bold')

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
end
end