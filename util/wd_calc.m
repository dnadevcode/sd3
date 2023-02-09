function [filt_log, wd] = wd_calc(sigma)

if nargin == 0
	sigma = 1;
	fprintf('No PSF given, assuming std dev of blurring = 1\n');
end
 

sixsig = round(6*sigma);
fsize = sixsig + (1-mod(sixsig,2)); % filter size - 6sigma rule
filt_log = fspecial('log', fsize, sigma); % log filter

if nargout > 1
    x0 = zeros(1,round(10*sigma));
    x1 = ones(1,round(10*sigma));
    x = [x0 x1];
    filt_gauss = fspecial('gauss',fsize,sigma); % gaussian filter
    x_blur = imfilter(x,filt_gauss,'replicate'); % signal blurred by gaussian
    x_log = imfilter(x_blur,filt_log,'replicate'); % applied laplacian of gaussian on the edge
    % xmid = x_log(round(size(x_log,1)/2),:); % middle column
    [~,xmin]= min(x_log);
    [~,xmax]= max(x_log);
    wd = abs(xmax-xmin); % distance between minimum and maximum of log filter
end
 
end