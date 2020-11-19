function [ k,b ] = get_line_parameters( bw )

% gets line parameters
[row, col] = ind2sub(size(bw), find(~isnan(bw)));
vals = [row col];

[coeff,~,~,~,~,mu]= pca(vals);
% Calculate the slope of the line
%	(since the y axis is reversed in image analysis, a minus sign is introduced
%  in order for the number k to represent what is normally thought of as the slope)
k = -coeff(1,1)/coeff(2,1);

% intercept of first column (1 is top row)
b = mu(1) + (mu(2) - 1) * k;

end

