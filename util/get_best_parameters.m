function [ maxcoef,pos,or ] = get_best_parameters( xcorrs)
% get_best_parameters

% input xcorrs
% output rezMax

% f gives the max over the two rows, s gives the indice of
% which row has the max (1 or 2)
[f,s] =max(xcorrs);

% sort the max scores, ix stores the original indices
[ b, ix ] = sort( f(:), 'descend' );

% choose the best three scores. (change this in the future?)
indx = b(1:3)' ;

% save the best three max score, and their orientation
maxcoef = indx;
or = s(ix(1:3)');
pos = ix(1:3)';

end

