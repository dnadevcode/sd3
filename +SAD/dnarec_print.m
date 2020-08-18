function printName = dnarec_print(output,experiment,actions,optics,runNo,sets)
nImage = numel(output);

% Extract all barcode lengths and dot numbers
imDots = zeros(1,nImage);
imBarLength = zeros(1,nImage);
imBars = zeros(1,nImage);
hasDots = zeros(1,nImage);
for i = 1:nImage
	hasDots(i) = isfield(output{i},'dots');
	nBars = numel(output{i}.expBars);
	imBars(i) = nBars;
	nDots = zeros(1,nBars);
	barLength = zeros(1,nBars);
	for j = 1:nBars
		if hasDots(i)
			nDots(j) = output{i}.dots{j}.N;
			imDots(i) = imDots(i) + nDots(j);	
		end
		barLength(j) = numel(output{i}.expBars{j}.rawBarcode)*optics.pixelSize/1000;
		imBarLength(i) = imBarLength(i) + barLength(j);
	end
end
% Calculate the average dot per micrometer in each image
imDotsPerLength = imDots./imBarLength;

% Initiate printing - make file with corrects filename - make new file if old is present
printName = print_version(nImage,experiment,output{1}.name,runNo);

% Print overall results
fid = fopen(printName,'w');
fprintf(fid,'Results for the analysis of %s\n',printName(1:end-12));
fprintf(fid,'\n Total number of barcodes: %i \n',sum(imBars));
fprintf(fid,'\n Total length of barcodes: %.1f micrometer \n',sum(imBarLength));
if hasDots
	fprintf(fid,'\n Total number of dots    : %i \n',sum(imDots));
	fprintf(fid,'\n Average dots/micron     : %.6f \n',sum(imDots)/sum(imBarLength));
end
fprintf(fid,'   Note that these number relate to the OBSERVED length of the molecule \n');
fprintf(fid,'   and might not accurately represent contour length! \n');
fprintf(fid,'----------------------------------------------------------------------- \n');
fprintf(fid,'----------------------------------------------------------------------- \n');
if nImage > 1
   % Print results for each image in folder
   for i = 1:nImage
      fprintf(fid,'\nFor image %s',output{i}.name);
      if actions.autoThreshBars
          molThresh = [num2str(output{i}.molScoreLim) '(Auto)'];
      else
          molThresh = [num2str(output{i}.molScoreLim) '(Manual)'];
      end
      fprintf(fid,'\n Minimum molecule score  :%s\n',molThresh); 
      fprintf(fid,'\n Total number of barcodes: %i \n',imBars(i));
      fprintf(fid,'\n Total length of barcodes: %.1f micrometer \n',imBarLength(i));
	  if hasDots(i)
          if actions.autoThreshDots
             dotThresh = [num2str(output{i}.dotScoreLim) '(Auto)'];
          else
             dotThresh = [num2str(output{i}.dotScoreLim) '(Manual)'];
          end
          fprintf(fid, '\n Minimum dot score      : %s \n',dotThresh);
	      fprintf(fid,'\n Total number of dots    : %i \n',imDots(i));
	      fprintf(fid,'\n Average dots/micron     : %.6f \n',imDotsPerLength(i));
	  end
	  fprintf(fid,'----------------------------------------------- \n');
   end
end
fprintf(fid,'Analysis settings:\n');
lengthLims = output{1}.lengthLims;
widthLims = output{1}.widthLims;
fprintf(fid,' Molecule length limits      : %.1f - %.1f pixels \n',lengthLims(1),lengthLims(2));
fprintf(fid,' Molecule width limits       : %.1f - %.1f pixels \n',widthLims(1),widthLims(2));
%fprintf(fid,'Molecule eccentricity limit : \n');
fprintf(fid,' Min. dot to end distance    : %.1f pixels \n',sets.dotMargin);
fprintf(fid,' Optics settings             : NA = %.2f, pixel size = %.2f, wavelength = %.2f. \n',optics.NA,optics.pixelSize,optics.waveLength);
% if actions.removeRegions
%     remSet = 'On';
% else
%     remSet = 'Off';
% end
% fprintf(fid,' Remove regions setting      : %s \n',remSet);
fclose(fid);


end

function printName = print_version(nImage,experiment,firstName,runNo)
	if nImage > 1
% 		slashes = strfind(experiment.targetFolder,'/');
% 		nameType = experiment.targetFolder(slashes(end-1)+1:slashes(end)-1);
        [~, nameType] = fileparts(experiment.targetFolder);
        folderName = experiment.targetFolder;
	else
		nameType = firstName;
        nameType = regexprep(nameType, '[.]\S{1,4}', '');
        folderName = fileparts(experiment.targetFolder);
	end
	nameType = [nameType,'results_run'];
%	content = dir(experiment.targetFolder);
%	numbers = zeros(1,numel(content));
%	for i = 1:numel(content)
%		matchName = strfind(content(i).name,nameType);
%		if ~isempty(matchName)
%			number = str2double(content(i).name(numel(nameType)+1));
%			isnum = ~isempty(number);
%			if isnum
%				numbers(i) = number;
%			end
%		end
%	end
%	version = max(numbers)+1;	
	version = runNo;	
	printName = fullfile(folderName, [nameType, num2str(version),'.txt']);
end



