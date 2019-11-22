function wd = wd_calc(sigma)

if nargin == 0
	sigma = 1;
	fprintf('No PSF given, assuming std dev of blurring = 1\n');
end

x0 = zeros(round(10*sigma));
x1 = ones(round(10*sigma));
x = [x0 x1];
sixsig = round(6*sigma);
fsize = sixsig + (1-mod(sixsig,2));
filt_gauss = fspecial('gauss',fsize,sigma);
filt_log = fspecial('log',fsize,sigma);
x_blur = imfilter(x,filt_gauss,'replicate');
x_log = imfilter(x_blur,filt_log,'replicate');
xmid = x_log(round(size(x_log,1)/2),:);
[~,xmin]= min(xmid);
[~,xmax]= max(xmid);
wd = abs(xmax-xmin);
