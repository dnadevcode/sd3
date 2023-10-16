function [name1,name2] = save_image_separate(noisyImage2D,noisyImageYoyo,nameext)
    % Save singleframe, supportable for optiscan

    % create format supportable by optiscan
    name1 = strcat(['output/' nameext 't_Scan_C=0.tif' ]);
    name2 = strcat(['output/' nameext 't_Scan_C=1.tif' ]);

    A = squeeze(noisyImageYoyo);
    B=squeeze(noisyImage2D);

        
    imwrite(uint16(round(A)),name1 );%'Resolution' ,[1,1]'Compression','none'
    imwrite(uint16(round(B)), name2 );

%     imwrite(uint16(round(A)),name)
%     imwrite(uint16(round(B)), name, 'writemode', 'append')


end

