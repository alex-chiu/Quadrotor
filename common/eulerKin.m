function [omegaB] = eulerKin(phi, theta, eF)
% eulerKin : Function that calculates necessary omegaB for a desired set of
%            Euler rates.
%
% INPUTS
%
% phi -------- x-direction Euler angle in radians
%
% theta ------ y-direction Euler angle in radians
%
% eF --------- 3-by-1 desired vector of Euler angle rates
%
% OUTPUTS
%
% omegaB ----- 3-by-1 vector of necessary body rotation
%
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 2/3/2022
%+==============================================================================+

S = 1 / cos(phi) .* [cos(phi)*cos(theta), 0, cos(phi)*sin(theta)
                    sin(phi)*sin(theta), cos(phi), -cos(theta)*sin(phi)
                    -sin(theta), 0, cos(theta)];

omegaB = S \ eF;

end