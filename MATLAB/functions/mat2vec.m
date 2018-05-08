%  MATLAB Function < mat2vec >
% 
%  Author:      Michele Facchinelli
%  Date:        5th May 2017
%  Purpose:     distribute columns of matrix to a variable number of outputs
%  Input:
%   - matrix:       matrix
%   - amount:       number of columns to distribute to each vector
%  Output:
%   - varargout:	vectors

function varargout = mat2vec(matrix,amount)

%...Give default value to amount
if nargin == 1
    if size(matrix,2) == 1
        matrix = matrix';
    end
    amount = ones(1,size(matrix,2));
end

%...Assign size to output
varargout = cell(1,size(amount,2));

%...Check if input is cell
cellinput = iscell(matrix);

%...Loop over columns
i = 0;
upper = cumsum(amount);
lower = upper-amount+1;
if ~cellinput
    for col = 1:size(amount,2)
        i = i+1;
        varargout{i} = matrix(:,lower(i):upper(i),:);
    end
else
    for col = 1:size(amount,2)
        i = i+1;
        varargout{i} = matrix{:,lower(i):upper(i),:};
    end
end