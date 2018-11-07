function weiner
D = '/home/ramya/Desktop/Projects/Mitosis_cancer/redoutput1';
out='/home/ramya/Desktop/Projects/Mitosis_cancer/weineroutput1';
S = dir(fullfile(D,'*.png')); 
for k = 1:numel(S)
    F = fullfile(D,S(k).name);
    I = imread(F);
    %imshow(I);
    [J,noise_out] = wiener2(I,[2 2])
    %noise_out
    imshow(J);
    baseFileName = sprintf('%s.png',S(k).name );
    fullFileName = fullfile(out, baseFileName); 
    imwrite(J, fullFileName);
    S(k).data = I; % optional, save data.
end
