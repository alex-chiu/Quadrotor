function [R] = matrixWrapper(col, n)
% matrixWrapper : Converts a column vector into an nxn matrix.
%
% INPUTS
%
% col -------- Input column vector
%
% n ---------- Dimension of output matrix
%
% OUTPUTS
%
% R ---------- n-by-n rotation matrix
%
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 2/2/2022
%+==============================================================================+

R = zeros(n, n);

for i = 1:n
    R(:, i) = col((i - 1) * n + 1:i * n, 1);
end

end