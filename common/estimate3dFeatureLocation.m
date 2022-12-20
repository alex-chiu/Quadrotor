function [rXIHat, Px] = estimate3dFeatureLocation(M, P)
% estimate3dFeatureLocation : Estimate the 3D coordinates of a feature point
%                             seen by two or more cameras with known pose.
%
%
% INPUTS
%
% M ---------- Structure with the following elements:
%
%       rxArray = 1xN cell array of measured positions of the feature point
%                 projection on the camera's image plane, in pixels.
%                 rxArray{i} is the 2x1 vector of coordinates of the feature
%                 point as measured by the ith camera.  To ensure the
%                 estimation problem is observable, N must satisfy N >= 2 and
%                 at least two cameras must be non-colinear.
%
%      RCIArray = 1xN cell array of I-to-camera-frame attitude matrices.
%                 RCIArray{i} is the 3x3 attitude matrix corresponding to the
%                 measurement rxArray{i}.
%
%       rcArray = 1xN cell array of camera center positions.  rcArray{i} is
%                 the 3x1 position of the camera center corresponding to the
%                 measurement rxArray{i}, expressed in the I frame in meters.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the quad's
%                 sensors, as defined in sensorParamsScript.m
%
% OUTPUTS
%
%
% rXIHat -------- 3x1 estimated location of the feature point expressed in I
%                 in meters.
%
% Px ------------ 3x3 error covariance matrix for the estimate rxIHat.
%
%+------------------------------------------------------------------------------+
% References: Lecture Notes 
%
% Author: Alex Chiu
%
% Last Edited: 4/18/2022
%+==============================================================================+  

N = length(M.rxArray);
HMat = zeros(2*N, 4);
RMat = zeros(2*N, 2*N);

for i = 1:N
    % Construct transformation matrix P = K*[RCI,t]
    t = -M.RCIArray{i} * M.rcArray{i};
    Pk = P.sensorParams.K*[M.RCIArray{i} t];
    HMat(2*i-1:2*i, :) = ...
        [M.rxArray{i}(1) * P.sensorParams.pixelSize * Pk(3, :) - Pk(1, :); ...
         M.rxArray{i}(2) * P.sensorParams.pixelSize * Pk(3, :) - Pk(2, :)];
    RMat(2*i-1:2*i, 2*i-1:2*i) = P.sensorParams.pixelSize^2 * P.sensorParams.Rc;
end

Hr = HMat(:, 1:3);
z = -HMat(:, 4);
rXIHat = Hr' / RMat * Hr \ Hr' / RMat * z;
Px = inv(Hr' / RMat * Hr);