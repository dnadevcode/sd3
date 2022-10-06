function [] = czi_to_tif(cziFold,bfcovertPath)
    % https://docs.openmicroscopy.org/bio-formats/5.7.1/users/comlinetools/index.html
    
    if nargin < 2
    %     cziFold 
        bfcovertPath = '..\..\..\postdoc\CODES\bftools\bftools\bfconvert';
    end
    data = dir(fullfile(cziFold,'*.czi'));
    % data = dir('C:\Users\Lenovo\postdoc\DATA\DOTS\images\images\*.czi');
    % data = dir('C:\Users\Lenovo\postdoc\DATA\DOTS\test3\*.czi');

    for i=1:length(data)
        name =fullfile(data(i).folder,data(i).name);
        command = strcat([bfcovertPath ' '  name ' ' strcat(name,'.tif')]);
        [a,b] = system(command);
        ch1 = imread(strcat(name,'.tif'),1);
        imwrite(ch1,strcat(name,'_C=1.tif'));
        ch2 = imread(strcat(name,'.tif'),2);
        imwrite(ch2,strcat(name,'_C=0.tif'));
    end

end

