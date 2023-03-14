function [pmin,allscores] = autothresh_random_bars(bgPixels,sets,bgMedian)
    % calc autothresh based on random barcodes taken from background pixels

%     barLen = 10000;
    dotBars{1} = bgPixels(randperm(length(bgPixels)));

    import Core.DotDetection.log_dot_detection;
    [peaks,allscores] = log_dot_detection(dotBars,sets,bgMedian);

    pmin = mean(allscores)+sets.numSigmasAutoThresh*std(allscores);
end

