function [logBgPixels,pmin] = autothresh_mean_Std(bgPixels,bgMedian)


    logBgPixels = bgPixels-bgMedian;% ;


    pmin = 2*nanstd(logBgPixels(:)); % for correctness. should take same number as when averaging for extracting barcode

end

