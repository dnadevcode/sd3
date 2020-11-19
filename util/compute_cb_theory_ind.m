function [ barcodePxRes,bitmaskPxRes,barcodeBpRes,bitmaskBpRes ] = compute_cb_theory_ind(theorySeq,rawL,optics,sets,experiment)

seqL = length(theorySeq);

% binding probability
[prob] = cb_theory(theorySeq,experiment.concN,experiment.concY,sets.model.yoyoBindingConstant,sets.model.netropsinBindingConstant,1000,1);

% bp to pixel convertion ration
meanBpExt_pixels = rawL/seqL;

% Convolve with a Gaussian
psfBp = optics.sigma / meanBpExt_pixels;
ker = fftshift(images.internal.createGaussianKernel(psfBp, length(prob)))';
multF=conj(fft(ker))';
barcodeBpRes = ifft(fft(prob).*multF);

% get the bitmask in bp resolution
bitmaskBpRes = ones(1,length(barcodeBpRes));
untrustedBp = round(psfBp*sets.deltaCut);
% in case of linearity edge pixels are bitmasked
if  sets.isLinear == 1
  if untrustedBp > length(bitmaskBpRes)
    fprintf('Warning: Untrusted region is longer than molecule\n')
    bitmaskBpRes = zeros(1,length(bitmaskBpRes));
  end
  bitmaskBpRes(1:untrustedBp) = 0;
  bitmaskBpRes(end-untrustedBp+1:end) = 0;
end

% convert to px resolution
barcodePxRes = bp2px(barcodeBpRes, meanBpExt_pixels);
if length(barcodePxRes) ~= rawL
  fprintf('bp2px routine failed and gave barcode with length %i instead of %i. Rescaling barcode \n',length(barcodePxRes),rawL);
  v = linspace(1,length(barcodePxRes),rawL);
  barcodePxRes = interp1(barcodePxRes,v);
end
bitmaskPxRes = ones(1,rawL);
bitmaskPxRes([1:min(rawL,round(optics.sigma*sets.deltaCut)),max(1,rawL-round(optics.sigma*sets.deltaCut))+1:rawL]) = 0;
end

