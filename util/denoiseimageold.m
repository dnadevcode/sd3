function [imOut] = denoise_image( im,minN,step )

    if nargin < 2 
        minN = 100;
        step = 30;
    end
    % in this function we cut im into smaller segments, but later we can
    % piece them back together
    
    import PD.Core.Extraction.detrend_2d;

    [xMax,yMax] = size(im);
    imOut = zeros(size(im,1),size(im,2));
    for j=1:step
        m = minN+j*10;
        ii = 1;
        ims = cell(ceil(xMax/m), ceil(yMax/m));
        imdsd =  cell(ceil(xMax/m), ceil(yMax/m));
        for i=1:m:xMax
            jj = 1;
            for kk=1:m:yMax
                ims{ii,jj} = im(i:min(i+m-1,end),kk:min(kk+m-1,end));
                imdsd{ii,jj} =  detrend_2d(ims{ii,jj});
                jj = jj+1;            
            end
            ii = ii+1;

        end
        imOut = imOut +cell2mat(imdsd)./step;
    
    end

end

