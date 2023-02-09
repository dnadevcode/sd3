function [] = czi_to_tif(cziFold,bfcovertPath)
    % https://docs.openmicroscopy.org/bio-formats/5.7.1/users/comlinetools/index.html
    
    if nargin < 2
    %     cziFold 
        bfcovertPath = 'bfconvert';
    end
  
    
    data = dir(fullfile(cziFold,'*.czi'));

    for i=1:length(data)
        name =fullfile(data(i).folder,data(i).name);
        
        [fd,fm,fe] = fileparts(name);
        nameNew = strrep(name,fe,'.tif');

        command = strcat([bfcovertPath ' '  name ' ' nameNew]);
        
        if exist(nameNew,'file')
            delete(nameNew); % in case already exists tif, remove 
        end
    
        [a,b] = system(command);
        ch1 = imread(nameNew,1);
        imwrite(ch1,strrep(name,fe,'_C=0.tif'));
        ch2 = imread(nameNew,2);
        imwrite(ch2,strrep(name,fe,'_C=1.tif'));
        
        delete(nameNew); % in case already exists tif, remove 

    end

end

