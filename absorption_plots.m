%Absorption plots

path_lengths = [0.5 0.63 1 1.5];

for i = 1:length(path_lengths)
    attenuation(:,i) = exp(-absorption * path_lengths(i));
    plot(wavelengths, attenuation(:,i));
    hold on
end
