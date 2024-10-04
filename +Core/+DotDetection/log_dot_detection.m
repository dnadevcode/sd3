function [peaks,allscores] = log_dot_detection(dotBars,sets,bgMedian)
    % LOG based dot detection

    %   Returns:
    %       peaks -
    %       allscores


  oldsig = sets.sigma;
  newsig = 669/509 * sets.sigma;   % TODO: Ask user to provide wavelength for second channel
%   newsig = sets.sigma/( 669/509);   % TODO: Ask user to provide wavelength for second channel

  w = 1;
  dots = cell(1, numel(dotBars));
  % Use optics to define LoG filter
  n = ceil(6 * newsig);
  n = n + 1 - mod(n, 2);
  range = (n - 1) / 2;
  idx = -range:range;
  filt = -1 / (pi * newsig^4) * (1 - idx.^2 / (2 * newsig^2)) .* exp(-idx.^2 / (2 * newsig^2)); 
  peaks = cell(1, numel(dotBars));
  allscores = [];
%   gaus = @(x,mu,sig,amp)amp*exp(-(((x-mu).^2)/(2*sig.^2)));

  % sumI = sum(cellfun(@(x) nansum(x), dotBars));
  % L = sum(cellfun(@(x) sum(numel(x)), dotBars));
  % averageI = sumI / L;
  % Filter images one by one and calculate a score
  for i = 1:numel(dotBars)

    if sum(~isnan(dotBars{i})) > sets.lengthLims(1)
      barLength = numel(dotBars{i});
      [initnan, endnan] = nanfind(dotBars{i});
      dotBars{i}(isnan(dotBars{i})) = 0; % keep only non-nans. Here can't be any nan's in the middle of barcode!
      dotBars{i} = dotBars{i}(1 + initnan:end - endnan)-bgMedian;
      dotBars{i}( dotBars{i} <0) = 0;
      %dotBars{i} = dotBars{i}(~isnan(dotBars{i}));
      logBar = imfilter(dotBars{i}, filt);
      [~, peaklocs] = findpeaks(-logBar); % find peaks after fuker

      %
      peaks{i}.scores = zeros(1, numel(peaklocs));
      peaks{i}.locations = peaklocs; % locations: shifted by initnan
      %peaks{i}.depth = min(peaklocs-initnan-oldsig,barLength-peaklocs-endnan-oldsig);
      peaks{i}.depth = max(zeros(1,length(peaklocs)),min(peaklocs - oldsig, barLength - peaklocs - endnan - initnan - oldsig)); % how far along barcode
      peaks{i}.leftOffset = initnan; % offsets
      peaks{i}.rightOffset = endnan;

      for j = 1:numel(peaklocs)
        idx = max(1, round(peaklocs(j) - w * newsig)):min(numel(dotBars{i}), round(peaklocs(j) + w * newsig));
        peaks{i}.scores(j) = nanmean(dotBars{i}(idx));% / bgMedian; % why
        peaks{i}.dotintensity(j) = sum(dotBars{i}(idx))/0.6827; % total dot intensity
%         divide by background? better substract?
% 
%         idx = @(x) max(1, round(peaklocs(j) - w * x)):min(numel(dotBars{i}), round(peaklocs(j) + w * x));
%         y = @(x) gaus(idx(x),peaklocs(j),x,dotBars{i}(peaklocs(j)));
%         minfun = @(x) sum(dotBars{i}(idx(x)))/sum(y(x))-1;
%         peaks{i}.bestSigma(j) = fzero(@(x) minfun(x),[0.1 5]);
      end

      allscores = [allscores peaks{i}.scores];
    else
      peaks{i}.scores = [];
      peaks{i}.dotintensity = [];
      peaks{i}.locations = [];
      peaks{i}.depth = [];
      peaks{i}.leftOffset = [];
      peaks{i}.rightOffset = [];
    end

  end

end

%%
% gaus = @(x,mu,sig,amp)amp*exp(-(((x-mu).^2)/(2*sig.^2)));
% idx = max(1, round(peaklocs(j) - w * newsig)):min(numel(dotBars{i}), round(peaklocs(j) + w * newsig));
% y = gaus(idx,peaklocs(j),newsig,dotBars{i}(peaklocs(j)));
% scaleFactor = sum(dotBars{i}(idx))/sum(y)
% 
% fmin = @(x) sum(dotBars{i}(idx))/sum(y)
% 
% sum(dotBars{i}(idx))/0.6827*scaleFactor
% 
% figure,plot(dotBars{i})
% hold on
% plot(idx,y)

% idx = @(newsig) max(1, round(peaklocs(j) - w * newsig)):min(numel(dotBars{i}), round(peaklocs(j) + w * newsig));
% y = @(newsig) gaus(idx(newsig),peaklocs(j),newsig,dotBars{i}(peaklocs(j)));
% minfun = @(x) sum(dotBars{i}(idx(x)))/sum(y(x))-1;
% bestSigma = fzero(@(x) minfun(x),1.5);
% 
% figure,plot(dotBars{i})
% hold on
% plot(idx(bestSigma),y(bestSigma))
% 
% sum(dotBars{i}(idx(bestSigma)))/0.6827

%%

function [initnan, endnan] = nanfind(bar)
  stillNaN = isnan(bar(1));
  i = 0;

  while stillNaN && i < numel(bar)
    i = i + 1;
    stillNaN = isnan(bar(i));
  end

  initnan = i;

  i = numel(bar);
  stillNaN = isnan(bar(i));

  while stillNaN && i > 1
    i = i - 1;
    stillNaN = isnan(bar(i));
  end

  endnan = numel(bar) - i;
end
