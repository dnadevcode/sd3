function [ randomSequences ] = generate_random_sequences(lengthOfSequences,optics,sets)
% generate_random_sequences
% Generates random sequences
%
%     Args:
%
%     Returns:
%         sets: Output settings

% first load pregenerated zero model
meanFFT = load(strcat([sets.pval.nulModelPath sets.pval.nulModelName]));
numberOfSequences = sets.pval.numRandBarcodes;

randomSequences = cell(1,numberOfSequences);

oldLength = length(meanFFT.meanFFTEst);

fftPr = meanFFT.meanFFTEst(1:floor((end+3)/2));

f1 = [0:length(1:(lengthOfSequences/2))]*(1/lengthOfSequences);

f2 = [0:length(fftPr)-1]*(1/oldLength);

% and the Parseval's identity also needs to be satisfied
len1 = oldLength*(oldLength-1);
len2 = lengthOfSequences*(lengthOfSequences-1);

newB = zeros(1, lengthOfSequences);
newB(1) = meanFFT.meanFFTEst(1)*lengthOfSequences/oldLength;

newf = interp1(f2(2:end),fftPr(2:end),f1(2:end),'linear','extrap');
%figure,plot(newf)

newB(2:length(newf)+1) = newf(1:end);
newB(lengthOfSequences-(1:length(newf))+1) =  newf(1:end);

%figure, plot(newB)
% the sums a and b have to be the same
a = sum(meanFFT.meanFFTEst.^2)/(len1); b = sum(newB.^2)/(len2);

% if they are not the same, we renormalise intVal values 2:end by
% constant konst
konst=sqrt((len2*a-newB(1).^2)/(len2*b-newB(1).^2));

% and define intValNorm
interpData = abs([newB(1) newB(2:end).*konst]);

% Why was once gaussian filter shifted and this one not?
ker = circshift(images.internal.createGaussianKernel(optics.sigma, lengthOfSequences),round(lengthOfSequences/2))';


halfL = floor(lengthOfSequences/2);

s = rng;
rng(0,'twister'); % generate reproducible random sequences
for i=1:numberOfSequences
  PR1 = exp(2i*pi*rand(1,halfL));
  PR2 = fliplr(conj(PR1));
  
  if mod(lengthOfSequences,2)==0
    PR = [1 PR1(1:end-1) 1 PR2(2:end)];
  else
    PR = [1 PR1 PR2];
  end
  
  randomSequences{i} = ifft((interpData.*PR).*conj(fft((ker))));
end
rng(s)

end

