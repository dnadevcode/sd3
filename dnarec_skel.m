function output = dnarec_skel(target)

% Store data path and assumptions on experimental conditions

if nargin == 0
 	experiment.targetFolder = '/Users/krog/Dropbox/dnabarcoding/SDD_dots/mytest/'; % Folder containing images and optics.txt file
else
	experiment.targetFolder = target;
end
targetFile = ''; % Optional specified target image. If the whole folder should be analysed, leave "targetFile = '';"
experiment.barFlag = 'C=1'; % Flag for specifying that the image is an image of barcodes. Only images containing the flag in their name will be analysed
experiment.dotFlag = 'C=90'; % Flag for specifying that the image is a dot image. The flag must be contanied in the image name.
experiment.opticsFile = 'optics.txt'; % Optionally specified optics file location (full path needed). If this is specified, the software assumes the optical setup to be identical for all images

sets.dotMargin = 4; % Minimum distance (in pixels) from molecule end for a dot to be included in the analysis.

%%% IFF only one (DNA molecule) image should be analyzed, write the full filename here
%targetFile = 'Snap-1283.czi - C=1.tif';
%%% else, comment the line.
	
% Specify functions to be used for analysis. Set these equal 1 to have the action performed or 0 to skip the action.
% You do not need to worry about these unless you would like to avoid showing or saving images / barcodes.

 actions.denoiseImage = 0; 			% Denoise the image using trend removal in 2D (shouldnt be needed for benchtop microscopy)
 actions.removeRegions = 0;			% Remove regions outside a central brigth region (can remove too much if illumination is uniform)
 actions.checkSameOrientation = 1;	% Check that images from the same stack (of the same subject) are identically oriented. Only applies or image stacks
 actions.makeSameSize = 1;			% Makes images in imagestacks the same size
 actions.showScores = 1; 			% Show a histogram of the scores for the regions as well as the dots, which helps setting tweak parameters
 actions.showMolecules = 1;			% Show plots of detected molecules and images
 actions.saveMolecules = 0;			% Save individual molecule and dot images (2D)
 actions.saveBars = 0; 				% Save the generated barcodes and dot barcodes
 actions.autoThreshBars = 0;
 actions.autoThreshDots = 0;
% actions.computeFreeConcentrations = 1;
 actions.getConsensus = 1;			% 


%%%% Note for users %%%%
% The rest of the parameters here are only used for automatically running
% p-value calculation and the consensus routine. They can be ignored if only generating barcode .tifs and movies of the corresponding molecules

 experiment.concY = 4e-2; 		% YOYO1 concentration in units of 1E-6 molar
 experiment.concN = 6; % Netropsin concentration in units of 1E-6 molar
% exp.concDNA = 0.2; 		% DNA concentration in units of 1E6 molar
% exp.seqence_path = PATH_OF_TARGET_SEQUENCE % Specify this if a custom sequences is to be analysed separately


% Specify settings (untrusted regions, model for theoretical barcodes etc.)
 sets.deltaCut = 1; % 
 import SAD.cb_model;
 sets.model = cb_model();
 sets.pval.nulModelPath = 'util/';
 sets.pval.nulModelName = 'meanFFT.mat';
 sets.pval.numRandBarcodes = 100;
% import SAD.import_sequence;
% sets.dnaSequence = import_sequence('./');
 sets.isLinear = 1;
 if isempty(targetFile)
	 import SAD.import_images_simple
	 functions.img_import = @() import_images_simple(experiment,actions);
 else
	 import SAD.import_single_image
	 functions.img_import = @() import_single_image(targetFile,experiment,actions);
 end 
 import SAD.denoise_images;
 functions.img_denoise = @(images) denoise_images(images,actions);
 functions.img_segment = @(images,imageNames,imageNumber,saveName) segment_image(images,imageNames,imageNumber,saveName,experiment,actions);
 import SAD.extract_barcodes;
 functions.bc_extract = @(movies,optics,lengthLims,imageNumber,saveName) extract_barcodes(movies,optics,lengthLims,imageNumber,saveName,sets,experiment,actions); 
 import SAD.analyse_barcodes;
 functions.bc_analyse = @(barcodes,optics,lengthLims) analyse_barcodes(barcodes,optics,lengthLims,sets,experiment);
 import SAD.generate_consensus;
 functions.consensus_gen = @(barcodes) generate_consensus(barcodes);

import SAD.dna_process_folder
output = dna_process_folder(experiment,functions,actions,sets);

save([experiment.targetFolder,'dnarecoutput'],'output')
save([experiment.targetFolder,'dnarecoutput'],'actions','-append')
save([experiment.targetFolder,'dnarecoutput'],'experiment','-append')
