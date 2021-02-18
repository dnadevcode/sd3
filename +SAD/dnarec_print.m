function printName = dnarec_print(output, experiment, actions, optics, runNo, sets)
  nImage = numel(output);
  fragLengthRangeBp = sets.fragLengthRangeBp;

  % Extract all barcode lengths and dot numbers
  imBarLengthAll = cell(1, nImage);
  imDotsAll = cell(1, nImage);
  imDotsRejected = cell(1, nImage);
  imBars = zeros(1, nImage);
  imHasDots = zeros(1, nImage);

  for i = 1:nImage
    imHasDots(i) = isfield(output{i}, 'dots');
    nBars = numel(output{i}.expBars);
    imBars(i) = nBars;
    nDots = zeros(1, nBars);
    barLength = zeros(1, nBars);
    rejectedDotBars = false(1, nBars);

    for j = 1:nBars

      if imHasDots(i)
        nDots(j) = output{i}.dots{j}.N;
        rejectedDotBars(j) = output{i}.dots{j}.rejected;
      end

      barLength(j) = numel(output{i}.expBars{j}.rawBarcode) * optics.pixelSize / 1000;
    end

    imBarLengthAll{i} = barLength;
    imDotsRejected{i} = rejectedDotBars;
    imDotsAll{i} = nDots;
  end

  imBarLength = cellfun(@sum, imBarLengthAll);
  imDotBarLength = cellfun(@(a, b) sum(a(not(b))), imBarLengthAll, imDotsRejected);
  imDots = cellfun(@(a, b) sum(a(not(b))), imDotsAll, imDotsRejected);
  % Calculate the average dot per micrometer in each image
  imDotsPerLength = imDots ./ imDotBarLength;

  % Initiate printing - make file with corrects filename - make new file if old is present
  printName = print_version(nImage, experiment, output{1}.name, runNo);

  % Print overall results
  fid = fopen(printName, 'w');
  fprintf(fid, 'Results for the analysis of %s\n', printName(1:strfind(printName, '.') - 1));
  fprintf(fid, '\n Total number of barcodes: %i \n', sum(imBars));
  fprintf(fid, '\n Total length of barcodes: %.1f micrometer \n', sum(imBarLength));
  fprintf(fid, '\n Average length of barcodes: %.1f micrometer \n', mean(horzcat(imBarLengthAll{:})));
  fprintf(fid, '\n Number of barcodes of lengths 0-%.1f micrometer: %.0f \n', ...
    fragLengthRangeBp(1), get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
    0, fragLengthRangeBp(1)));

  if length(fragLengthRangeBp) > 1

    for j = 2:length(fragLengthRangeBp)
      fprintf(fid, ' Number of barcodes of lengths %.1f-%.1f micrometer: %.0f \n', ...
        fragLengthRangeBp(j - 1), fragLengthRangeBp(j), ...
        get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
        fragLengthRangeBp(j - 1), fragLengthRangeBp(j)));
    end

  end

  fprintf(fid, ' Number of barcodes of lengths >%.1f micrometer: %.0f \n', ...
    fragLengthRangeBp(end), get_num_frags_in_range(horzcat(imBarLengthAll{:}), ...
    fragLengthRangeBp(end), inf));

  if any(imHasDots)
    fprintf(fid, '\n Total number of dots    : %i \n', sum(imDots));
    fprintf(fid, '\n Average dots/micron     : %.6f \n', sum(imDots) / sum(imDotBarLength));
  end

  fprintf(fid, '   Note that these number relate to the OBSERVED length of the molecule \n');
  fprintf(fid, '   and might not accurately represent contour length! \n');
  fprintf(fid, '----------------------------------------------------------------------- \n');
  fprintf(fid, '----------------------------------------------------------------------- \n');

  if nImage > 1
    % Print results for each image in folder
    for i = 1:nImage
      fprintf(fid, '\nFor image %s', output{i}.name);

      if actions.autoThreshBars
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

        if actions.autoThreshDots
          dotThresh = [num2str(output{i}.dotScoreLim) '(Auto)'];
        else
          dotThresh = [num2str(output{i}.dotScoreLim) '(Manual)'];
        end

        fprintf(fid, '\n Minimum dot score       : %s \n', dotThresh);
        fprintf(fid, '\n Total number of dots    : %i \n', imDots(i));
        fprintf(fid, '\n Average dots/micron     : %.6f \n', imDotsPerLength(i));
      end

      fprintf(fid, '----------------------------------------------- \n');
    end

  end

  fprintf(fid, 'Analysis settings:\n');
  lengthLims = output{1}.lengthLims;
  widthLims = output{1}.widthLims;
  fprintf(fid, ' Minimum molecule score      : %.3g (User set value) \n', experiment.lowLim);
  fprintf(fid, ' Minimum dot score           : %.3g (User set value) \n', experiment.dotScoreMin);
  fprintf(fid, ' Molecule length limits      : %.1f - %.1f pixels \n', lengthLims(1), lengthLims(2));
  fprintf(fid, ' Molecule width limits       : %.1f - %.1f pixels \n', widthLims(1), widthLims(2));
  fprintf(fid, ' Molecule eccentricity limit : %.2f \n', experiment.elim);
  fprintf(fid, ' Min. mol-to-convex-hull     : %.2f \n', experiment.ratlim);
  fprintf(fid, ' Min. dot to end distance    : %.1f pixels \n', sets.dotMargin);
  fprintf(fid, [' Optics settings             : NA = %.2f,' ...
              ' pixel size = %.2f nm, \n' ...
              '   wavelength = %.2f nm, sigma_psf = %.2f nm, sigma_LoG = %.2f nm. \n'], ...
    optics.NA, optics.pixelSize, optics.waveLength, ...
    optics.sigma * optics.pixelSize, optics.logSigma * optics.pixelSize);
  fclose(fid);

end

function printName = print_version(nImage, experiment, firstName, runNo)

  if nImage > 1
    % [~, nameType] = fileparts(experiment.targetFolder);
    folderName = experiment.targetFolder;
  else
    % nameType = firstName;
    % nameType = regexprep(nameType, '[.]\S{1,4}', '');
    folderName = fileparts(experiment.targetFolder);
  end

  % nameType = [nameType,'results_run'];
  nameType = 'results_run';
  version = runNo;
  printName = fullfile(folderName, [nameType, num2str(version), '.txt']);
end

function numFrags = get_num_frags_in_range(barLengths, lower, upper)
  numFrags = length(barLengths(barLengths > lower & barLengths <= upper));
end
