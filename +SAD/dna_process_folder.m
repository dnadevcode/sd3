function output = dna_process_folder(experiment,functions,actions,sets)

% This routine analyses the image contents of the folder specified in the "exp" struct
% via the functions specified in the "functions" struct according to te actions
% specified in the "actions" struct.

addpath('util/');
% addpath('util/lldev/src/MATLAB/');
import SAD.dnarec_print

% Perform preliminary check to see if all required files are accessible
% check = prel_check(experiment,functions,actions);
%
% if check
% Import and crop images
[images,imageNames] = functions.img_import();
output = cell(1,numel(images));

if prompt_figure_excess(length(images), sum(cellfun(@(x) isfield(x, 'dotIm'), images)), actions)
    return
end

% Specify save name
% slashes = strfind(experiment.targetFolder,'/');
% folderName = experiment.targetFolder(slashes(end-1)+1:slashes(end)-1);
folderName = subsref(dir(experiment.targetFolder), substruct('.', 'folder'));
runNo = 0;
outputExist = 1;
while outputExist
    runNo = runNo + 1;
%     if runNo > 98
%         fprintf('Only %i run outputs supported, overwriting last result (run %i).\n',runNo,runNo);
%         break;
%     end
    barName = fullfile(folderName, ['barcodes_run', num2str(runNo)]);
    dotName = fullfile(folderName, ['dotbars_run', num2str(runNo)]);
    molName = fullfile(folderName, ['molecules_run', num2str(runNo)]);
    resName = fullfile(folderName, ['results_run', num2str(runNo), '.txt']);
    barcodeExist = isfolder(barName);
    dotbarExist = isfolder(dotName);
    moleculesExist = isfolder(molName);
    resExist = isfolder(resName);
    outputExist = barcodeExist || dotbarExist || moleculesExist || resExist;
end
if actions.saveMolecules
    mkdir(molName)
end
if actions.saveBars
    mkdir(barName)
    mkdir(dotName)
end
for i = 1:numel(images)
    fprintf('\nAnalysing image %s.\n',imageNames{i});
    % Denoise images and remove artifacts
    cleanImages = functions.img_denoise(images{i}.registeredIm);
    if isfield(images{i},'dotIm')
        cleanImages.dotIm = images{i}.dotIm;
    end
    
    % Segment image
    [movies,~,optics,lengthLims,molScoreLim,widthLims] = functions.img_segment(cleanImages,imageNames{i},i,runNo);
    
    % Extract barcodes
    [barcodes,dotScoreLim] = functions.bc_extract(movies,optics,lengthLims,i,runNo);
    
    if actions.showMolecules
        % Mark barcodes in molecule image
        mark_bars(movies,barcodes);
    end
    
    % Calculate p-values for specific sequence
    %	output = functions.bc_analyse(barcodes,optics);i
    output{i} = barcodes;
    output{i}.name = imageNames{i};
    output{i}.lengthLims = lengthLims;
    output{i}.widthLims = widthLims;
    output{i}.molScoreLim = molScoreLim;
    output{i}.dotScoreLim = dotScoreLim;
    output{i}.median = median(cleanImages.imAverage(:));
end
import SAD.dnarec_print
resultsName = dnarec_print(output,experiment,actions,optics,runNo,sets);
fprintf('\n-------------------------------------------------------------------\n');
fprintf('Analysis complete\n');
fprintf('Results saved in %s',resultsName);
fprintf('\n-------------------------------------------------------------------\n');

% Run the consensus barcode routine on the barcodes
%	consensusStruct = functions.consensus_gen(barcodes);

% Print results
%	dnarec_print(output)
% Mark significants
%	mark_matches(output,movies.pos,barcodes.delid,0,0.09)

%	% Save results
%	dnarec_save(images,cleanImages,movies,scores,barcodes,output,experiment,consensusStruct);
% else
% 	fprintf('Preliminary check failed. Analysis aborted.\n');
% 	output = NaN;
end
