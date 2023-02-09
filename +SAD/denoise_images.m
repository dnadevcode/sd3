function cleanImages = denoise_images(registeredIm,sets)

tic;
% if actions.denoiseImage == 1
% 	% Denoisify
%  	[imDenoised, imAverage ]= get_denoised_image(registeredIm);
%  else
 	imAverage = averageImages(registeredIm);
% imAverage = double(registeredIm);
% imDenoised = imAverage; % This is NOT clean
%  end

%% Remove circles
optics.logSigma = sets.logSigmaNm / sets.pxnm;
% dilateSize = 2;
dilateSize = 0;
rLims = [ceil(6 * optics.logSigma) 23];
  
resizefactor = 1/10; % resize factor in order to speed up calculation
rLims(1) = round(rLims(1)*resizefactor)+dilateSize;

% original size
sz = round(size(imAverage)*resizefactor);
% scale image down:
imScaled = imresize(imAverage,sz);

% expect the circles to have higher intensity than all else
lim = multithresh(imScaled,2);

circles = imScaled > lim(2);

% se = strel('disk',dilateSize); % dilate circles by a disk
% dilatedBW = imdilate(circles,se);

[centers, ~] = imfindcircles(circles, rLims);

% ww = (not(circles)).*(imAverage);
%figure,imshow(ww,[])

% nw = ww;
% nw(nw==0) = nan;


% [~,centers] = get_remove_circles(imDenoised, sets);

cleanImages.registeredIm = registeredIm;
cleanImages.imAverage = imAverage;
% cleanImages.imDenoised = imDenoised;
cleanImages.centers = round(centers./resizefactor);

% 
% figure,imagesc(imAverage)
% hold on
% for i=1:size(cleanImages.centers,1)
%     plot(cleanImages.centers(i,1),cleanImages.centers(i,2),'redx')
% end


t=toc;
fprintf('Detected %.1d circles. \n',size(cleanImages.centers,1));
fprintf('Image denoising completed in %.1f seconds. \n',t);
