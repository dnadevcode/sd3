function cleanImages = denoise_images(registeredIm,actions)

tic;
if actions.denoiseImage == 1
	% Denoisify
 	[imDenoised, imAverage ]= get_denoised_image(registeredIm);
 else
 	imAverage = averageImages(registeredIm);
 	imDenoised = imAverage; % This is NOT clean
 end

% Remove circles
[~,centers] = get_remove_circles(imDenoised, imAverage);

cleanImages.registeredIm = registeredIm;
cleanImages.imAverage = imAverage;
cleanImages.imDenoised = imDenoised;
cleanImages.centers = centers;
t=toc;
fprintf('Image denoising completed in %.1f seconds. \n',t);
