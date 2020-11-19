function [ parameters ] = compute_evd_params( cc,x0 )
%compute_evd_params
% Compute functional fomr parameters for max cc histogram
%     Args:
%       cc,x0
%     Returns:
%         parameters

f = @(y) abs(functional_form_helping_fun(y,cc));

% first parameter, the value should be found by running minimization
% in the interval [2,x0]
x2 = fminbnd(f,2,x0); % This is \nu

% second parameter
m = length(cc);
cc = min(cc,1);

denom = (-1/m*sum(log(1+betainc(cc.^2,1/2,x2/2-1)))+log(2));
N2 = 1/denom; % This is \lambda

parameters = [x2 N2];
end

