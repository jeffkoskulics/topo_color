function [rms, std_dev, average] = check_residual(guess,reference)
%check_residual outputs difference measurements of the residual either
%between subsequent iterations or between the iterations and the "real"
%data

diff = reference - guess;

if size(diff,3) == 3 %Surface normals (or colors), z component ignored 
                     % (or blue channel) ignored
    %RMS
    rms(1:2) = sqrt(sum(sum(diff(:,:,1:2).^2))/(size(diff,1)*size(diff,2)));

    %Std_dev
    std_dev(1) = std(reshape(diff(:,:,1),size(diff,1)*size(diff,2),1));
    std_dev(2) = std(reshape(diff(:,:,2),size(diff,1)*size(diff,2),1));
    
    %Mean
    average(1:2) = mean(mean(diff(:,:,1:2)));


else %Elevations

    %RMS
    rms = sqrt(sum(sum(diff.^2)))/(size(diff,1)*size(diff,2));

    %Std_dev
    std_dev = std(reshape(diff,size(diff,1)*size(diff,2),1));
    
    %Mean
    average = mean(mean(diff));
end