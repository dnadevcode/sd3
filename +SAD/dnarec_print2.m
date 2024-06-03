function printName = dnarec_print2(output, sets, runNo,filtered)
    % prints results as txt
    
    if nargin < 7
        filtered = 0;
    end
  nImage = numel(output);
  fragLengthRangeBp = sets.fragLengthRangeBp;

  % Extract all barcode lengths and dot numbers and vals
  imBarLengthAll = cell(1, nImage);
  imDotsAll = cell(1, nImage);
  imDotsRejected = cell(1, nImage);
  imDotsValsAll = cell(1, nImage);

  imBars = zeros(1, nImage);
  imHasDots = zeros(1, nImage);

  for i = 1:nImage
    imHasDots(i) = isfield(output{i}, 'dots2');
    nBars = numel(output{i}.expBars);
    imBars(i) = nBars;
    nDots = zeros(1, nBars);
    barLength = zeros(1, nBars);
    rejectedDotBars = false(1, nBars);
    dotInt = zeros(1, nBars);

    for j = 1:nBars

      if imHasDots(i)
        nDots(j) = output{i}.dots2{j}.N;
        rejectedDotBars(j) = output{i}.dots2{j}.rejected;
        dotInt(j) = sum(output{i}.dots2{j}.val);
      end

      barLength(j) = numel(output{i}.expBars{j}.rawBarcode) * sets.pixelSize / 1000;
    end

    imDotsValsAll{i} = dotInt;
    imBarLengthAll{i} = barLength;
    imDotsRejected{i} = rejectedDotBars;
    imDotsAll{i} = nDots;
  end
  
    dotIntAll  = cellfun(@sum, imDotsValsAll);
    
    imBarLength = cellfun(@sum, imBarLengthAll);
    imDotBarLength = cellfun(@(a, b) sum(a(not(b))), imBarLengthAll, imDotsRejected);
    imDots = cellfun(@(a, b) sum(a(not(b))), imDotsAll, imDotsRejected);
    allDots = horzcat(imDotsAll{:});
    allDotsRej = horzcat(imDotsRejected{:});

    allBarlength = horzcat(imBarLengthAll{:}); 
    %     imDotsBar = cellfun(@(a, b) a(not(b)), imDotsAll, imDotsRejected);
    
    % Calculate the average dot per micrometer in each image
    imDotsPerLength = imDots ./ imDotBarLength;
    
    % Initiate printing - make file with corrects filename - make new file if old is present
    printName = print_version(nImage, sets, output{1}.name, runNo, filtered);


  % Print overall results
  fid = fopen(printName, 'w');
  fprintf(fid, 'Results for the analysis of %s\n', printName(1:strfind(printName, '.') - 1));
  fprintf(fid, '\n Total number of barcodes: %i \n', sum(imBars));
  fprintf(fid, '\n Total length of barcodes: %.1f micrometer \n', sum(imBarLength));
  fprintf(fid, '\n Average length of barcodes: %.1f micrometer \n', mean(horzcat(imBarLengthAll{:})));
    [numBars,locations] = get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
    0, fragLengthRangeBp(1));
  fprintf(fid, '\n Number of barcodes of lengths 0-%.1f micrometer: %.0f ( %.3f dots/micron)\n', ...
    fragLengthRangeBp(1), numBars, sum(allDots(logical(locations.*~(allDotsRej))))/sum(allBarlength(locations)));

  if length(fragLengthRangeBp) > 1

    for j = 2:length(fragLengthRangeBp)
        [numBars,locations] = get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
        fragLengthRangeBp(j - 1), fragLengthRangeBp(j));
        fprintf(fid, ' Number of barcodes of lengths %.1f-%.1f micrometer: %.0f ( %.3f dots/micron) \n', ...
        fragLengthRangeBp(j - 1), fragLengthRangeBp(j), ...
        numBars,...
        sum(allDots(logical(locations.*~(allDotsRej))))/sum(allBarlength(locations))); % dots/micron
    end

  end

    [numBars,locations] = get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
    fragLengthRangeBp(end), inf);
    fprintf(fid, ' Number of barcodes of lengths >%.1f micrometer: %.0f ( %.3f dots/micron) \n', ...
    fragLengthRangeBp(end), numBars, sum(allDots(logical(locations.*~(allDotsRej))))/sum(allBarlength(locations)));

  if any(imHasDots)
    fprintf(fid, '\n Total number of dots    : %i \n', sum(imDots));
    fprintf(fid, '\n Total intensity of dots   : %i \n', sum(dotIntAll));

    fprintf(fid, '\n Average dots/micron     : %.6f \n', sum(imDots) / sum(imDotBarLength));
  end

  fprintf(fid, '   Note that these number relate to the OBSERVED length of the molecule \n');
  fprintf(fid, '   and might not accurately represent contour length! \n');
  fprintf(fid, '----------------------------------------------------------------------- \n');
  fprintf(fid, '----------------------------------------------------------------------- \n');

  if nImage >= 1
    % Print results for each image in folder
    for i = 1:nImage
      fprintf(fid, '\nFor image %s', output{i}.name);

      if sets.autoThreshBars
        molThresh = [num2str(output{i}.molScoreLim) '(Auto)'];
      else
        molThresh = [num2str(output{i}.molScoreLim) '(Manual)'];
      end

      fprintf(fid, '\n Minimum molecule score  : %s \n', molThresh);
      fprintf(fid, '\n Total number of barcodes: %i \n', imBars(i));
      fprintf(fid, '\n Total length of barcodes: %.1f micrometer \n', imBarLength(i));
      fprintf(fid, '\n Average length of barcodes: %.1f micrometer \n', mean(imBarLengthAll{i}));
      fprintf(fid, '\n Number of barcodes of lengths 0-%.1f micrometer: %.0f \n', ...
        fragLengthRangeBp(1), get_num_frags_in_range(horzcat(imBarLengthAll{i}), ...
        0, fragLengthRangeBp(1)));

      if length(fragLengthRangeBp) > 1

        for j = 2:length(fragLengthRangeBp)
          fprintf(fid, ' Number of barcodes of lengths %.1f-%.1f micrometer: %.0f \n', ...
            fragLengthRangeBp(j - 1), fragLengthRangeBp(j), ...
            get_num_frags_in_range(horzcat(imBarLengthAll{i}), ...
            fragLengthRangeBp(j - 1), fragLengthRangeBp(j)));
        end

      end

      fprintf(fid, ' Number of barcodes of lengths >%.1f micrometer: %.0f \n', ...
        fragLengthRangeBp(end), get_num_frags_in_range(horzcat(imBarLengthAll{i}), ...
        fragLengthRangeBp(end), inf));

      if imHasDots(i)

        if sets.autoThreshDots
          dotThresh = [num2str(output{i}.dotScoreLim2) '(Auto)'];
        else
          dotThresh = [num2str(output{i}.dotScoreLim2) '(Manual)'];
        end

        fprintf(fid, '\n Minimum dot score       : %s \n', dotThresh);
        fprintf(fid, '\n Total number of dots    : %i \n', imDots(i));
        fprintf(fid, '\n Average dots/micron     : %.6f \n', imDotsPerLength(i));
        fprintf(fid, '\n Intensity of dots    : %.6f \n', dotIntAll(i));
      
    
        for jj=1:length(output{i}.dots2)
            fprintf(fid, strcat(['\n Mol ' num2str(jj) ' [micrometer length] : %4.3f']),imBarLengthAll{i}(jj));
            fprintf(fid, strcat(['\n Number dots : %4d']),length(output{i}.dots2{jj}.locations));
            try
                fprintf(fid, strcat(['\n Dot score : %s ']),regexprep(num2str(output{i}.dots2{jj}.val),'\s+',','));
                fprintf(fid, strcat(['\n Depth [micrometer] : %s ']),regexprep(num2str(output{i}.dots2{jj}.depth*(optics.pixelSize / 1000)),'\s+',','));
            catch
            end
            fprintf(fid, '\n');

        end
      end
      fprintf(fid, '----------------------------------------------- \n');

%             fprintf(fid, '\n Total intensity of dots   : %i \n', sum(dotIntAll));


    end

  end

  % save all settings

  fprintf(fid, 'Analysis settings:\n');
  lengthLims = output{1}.lengthLims;
  widthLims = output{1}.widthLims;
  fprintf(fid, ' Pixel size    : %.2f nm (User set value) \n',  sets.pixelSize);
  fprintf(fid, ' Width of LoG filter (nm)    : %.2f nm (User set value) \n',  sets.logSigma * sets.pixelSize);
  %
  fprintf(fid, ' Molecule image flag    : %s (User set value) \n',  sets.barFlag);
  fprintf(fid, ' Dot image flag   : %s (User set value) \n',  sets.dotFlag);
  %
  fprintf(fid, ' Minimum log(EdgeScore)      : %.3g (User set value) \n', log(sets.lowLim));
  fprintf(fid, ' Minimum dot score           : %.3g (User set value) \n', sets.dotScoreMin);
  fprintf(fid, ' Molecule length limits      : %.1f - %.1f pixels \n', lengthLims(1), lengthLims(2));
  fprintf(fid, ' Molecule width limits       : %.1f - %.1f pixels \n', widthLims(1), widthLims(2));
  fprintf(fid, ' Molecule eccentricity limit : %.2f \n', sets.elim);
  fprintf(fid, ' Min. mol-to-convex-hull     : %.2f \n', sets.ratlim);
  fprintf(fid, ' Min. dot to end distance    : %.1f pixels \n', sets.dotMargin);
  fprintf(fid, ' Mol extraction method    : %.1f (1-line, 2-spline)\n', sets.extractionMethod);
  fprintf(fid, ' autoThreshBars    : %.2f  (User set value) \n',  sets.autoThreshBars);
  fprintf(fid, ' autoThreshDots    : %.2f  (User set value) \n',  sets.autoThreshDots);

  fprintf(fid, ' showScores    : %.2f  (User set value) \n',  sets.showScores);
  fprintf(fid, ' showMolecules    : %.2f  (User set value) \n',  sets.showMolecules);
  fprintf(fid, ' saveMolecules    : %.2f  (User set value) \n',  sets.saveMolecules);
  fprintf(fid, ' saveBars    : %.2f  (User set value) \n',  sets.saveBars);
  fprintf(fid, ' remove non-uniform noise    : %.2f  (User set value) \n',  sets.denoiseImages);

  fprintf(fid, [' Optics settings             : NA = %.2f,' ...
              ' pixel size = %.2f nm, \n' ...
              '   wavelength = %.2f nm, sigma_psf = %.2f nm, sigma_LoG = %.2f nm. \n'], ...
    sets.NA, sets.pixelSize, sets.waveLength, ...
    sets.sigma * sets.pixelSize, sets.logSigma * sets.pixelSize);
  fclose(fid);

end

function printName = print_version(nImage, experiment, firstName, runNo,filtered)

  if nImage > 1
    % [~, nameType] = fileparts(experiment.targetFolder);
    folderName = experiment.targetFolder;
  else
    % nameType = firstName;
    % nameType = regexprep(nameType, '[.]\S{1,4}', '');
    folderName = fileparts(experiment.targetFolder);
  end
  
  if nargin < 5
      filtered = 0;
  end

  % nameType = [nameType,'results_run'];
  if filtered
        nameType = 'results_filtered';
  else
      nameType = 'results_run';
  end
  version = runNo;
  printName = fullfile(folderName, [nameType, num2str(version),'-2', '.txt']);
end

function [numFrags,locations] = get_num_frags_in_range(barLengths, lower, upper)
    locations = barLengths > lower & barLengths <= upper;
    numFrags = length(barLengths(locations));
end
