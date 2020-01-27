function printName = dnarec_print(output,experiment,optics,runNo)
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
      fprintf(fid,'\n Total number of barcodes: %i \n',imBars(i));
      fprintf(fid,'\n Total length of barcodes: %.1f micrometer \n',imBarLength(i));
	  if hasDots(i)
	      fprintf(fid,'\n Total number of dots    : %i \n',imDots(i));
	      fprintf(fid,'\n Average dots/micron     : %.6f \n',imDotsPerLength(i));
	  end
	  fprintf(fid,'----------------------------------------------- \n');
   end
end
fclose(fid);


end

function printName = print_version(nImage,experiment,firstName,runNo)
	if nImage > 1
		slashes = strfind(experiment.targetFolder,'/');
		nameType = experiment.targetFolder(slashes(end-1)+1:slashes(end)-1);
	else
		nameType = firstName;
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
	printName = [experiment.targetFolder,nameType,num2str(version),'.txt'];
end



