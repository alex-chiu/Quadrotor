function [col] = matrixUnwrapper(R)
% matrixWrapper : Converts a matrix into a column vector
%
% INPUTS
%
% R ---------- Input matrix
%
% n ---------- Dimension of output matrix
%
% OUTPUTS
%
% col -------- Output column vector
%
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 2/2/2022
%+==============================================================================+

n = length(R);
col = zeros(n, 1);

for i = 1:n
    col((i - 1) * n + 1:i * n, 1) = R(:, i);
end

end