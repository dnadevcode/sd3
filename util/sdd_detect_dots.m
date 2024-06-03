function [dots, pmin] = sdd_detect_dots(dotBars, sets, imageNumber, imageName, bgPixels,tiles)
    
    %   Args:
    %       dotBars

    %   Return:
    %       dots - dot structure
    %       pmin - 

    pmin = nan;
    bgMedian = median(bgPixels(:),'omitnan');

    import Core.DotDetection.log_dot_detection;
    sets.dotdetMethod = 'log';
    switch sets.dotdetMethod
        case 'log'
            [peaks, allscores] = log_dot_detection(dotBars,sets,bgMedian);
        case 'localmaxima' %todo
        case 'findpeaks'
        case 'sfw'
        otherwise
    end

    import Core.DotDetection.autothresh_mean_Std;
    import Core.DotDetection.autothresh_random_bars;

%     sets.autoThreshDotsMethod = 'meanstd';
    if sets.autoThreshDots % ?? 
        switch sets.autoThreshDotsMethod
            case 'meanstd'
                [logBgPixels,pmin] = autothresh_mean_Std(bgPixels,bgMedian);
            case 'randomBars'
                [pmin, logBgPixels] = autothresh_random_bars(bgPixels,sets,bgMedian);
            otherwise
        end
        if ~isempty(tiles)
            axes(tiles.bgDotScores);
        %     figure(5 + (imageNumber - 1) * 5)
            hbg = histogram(logBgPixels(:), 20);
            title([imageName, ' background (dot) score histogram'])
        end
    else
         pmin = sets.dotScoreMin;
  end

  medint = median(allscores(allscores > pmin));

  if sets.showScores
%     t = tiledlayout(hPanelResult,2,2,'TileSpacing','compact');
    axes(tiles.dotScores);
%     figure(5 + (imageNumber - 1) * 5)
    h2 = histogram(allscores, 20);
    title([imageName, ' dot scores'])
    hold on
    %    line([autoThresh autoThresh],[0 max(h2.Values)],'LineStyle','--','Color','red','LineWidth',2)
    if sets.autoThreshDots
      mess = 'Automated threshold';
      line([pmin pmin], [0 max(h2.Values)], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 2)
    else
      mess = 'User set threshold';
      line([pmin pmin], [0 max(h2.Values)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 2)
    end

    text(1.1 * pmin, 2/3 * max(h2.Values), mess, 'FontSize', 14)
    hold off
    
%     % Alternative: all intensities after filt
%     axes(tiles.dotScoresFilt);
%     %     figure(5 + (imageNumber - 1) * 5)
%     hfilt = histogram(allscores(allscores >pmin), 20);
%     title([imageName, ' dot scores filtered'])
  end

  % Locate peak positions in images with "high" scores
  totdots = 0;
  marginDots = 0;
  totInt = 0; % total intensity

  dots = cell(1,numel(dotBars));
  for i = 1:numel(dotBars)
    mask = peaks{i}.scores > pmin;
    endMask = peaks{i}.depth >= sets.dotMargin;
    accMask = and(mask, endMask);
    dots{i}.marginDots = sum(mask .*~endMask);
    marginDots = dots{i}.marginDots + marginDots;

    dots{i}.locations = peaks{i}.locations(accMask);
    dots{i}.depth = peaks{i}.depth(accMask);
    dots{i}.N = sum(accMask);
    totdots = totdots + dots{i}.N;
    dots{i}.val = peaks{i}.scores(accMask);% / medint;% ??
    dots{i}.leftOffset = peaks{i}.leftOffset;
    dots{i}.rightOffset = peaks{i}.rightOffset;

    dots{i}.rejected = length(dotBars{i}) <= 2*sets.dotMargin;
    totInt = totInt + sum(dots{i}.val);

    if dots{i}.N > 0 && sets.showDotPeaks
      logBar = -imfilter(dotBars{i}, filt);
      thisScores = peaks{i}.scores(accMask);
%       axes(tiles.dotPeaks);
        figure
      plot(dotBars{i});
      hold on
      plot(logBar / mean(logBar) * mean(dotBars{i}));

      for j = 1:length(dots{i}.locations)
        loc = dots{i}.locations(j);
        xline(loc, '--r');
        text(loc, dotBars{i}(loc), num2str(thisScores(j), 2));
      end

      %     scatter(dots{i}.locations, dotBars{i}(dots{i}.locations), 'ro');
      hold off
      xlabel('Position on molecule (pixels)')
      ylabel('Rescaled intensity')
      legend({'Dot-bar'; 'LoG-filtered Dot-bar'; 'Detected dots'}, ...
        'Location', 'northeastoutside')
    end

  end

  fprintf('Found %i dots with total intensity (after bg substr) %f.\n', totdots, totInt);
  fprintf('Rejected %i dots due to vicinity to molecule end.\n', marginDots);
end
