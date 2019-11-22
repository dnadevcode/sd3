function [ imDenoised,imAverage ] = get_denoised_image( registeredIm)

import PD.Core.Extraction.compare_simple_to_average;
% compare one frame vs average
[imAverage] = compare_simple_to_average( registeredIm );

import PD.Core.Extraction.denoise_image;
[imDenoised] = denoise_image( imAverage );

import PD.Core.Extraction.compare_2;
%compare_2(imDenoised,imAverage);
end

