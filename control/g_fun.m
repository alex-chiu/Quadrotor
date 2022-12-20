function [rCL,RCL] = g_fun(x,P)
% g_fun : Nonlinear function that maps the inputs to a full camera pose
%         expressed in the local ENU frame whose origin is P.sensorParams.r0G
%
% INPUTS
%
% x ---------- 6x1 vector arranged as [azRad;elRad;rotRad;rP], where
%
%     azRad -- Azimuth angle for the unit vector pointing from the
%              secondary to the primary antenna (i.e., toward the camera),
%              relative to the local ENU frame, in rad.  Azimuth is defined
%              as zero toward local north and increases clockwise.
%
%     elRad -- Elevation angle for the unit vector pointing from the
%              secondary to the primary antenna (i.e., toward the camera),
%              relative to the local ENU frame, in rad.  Elevation is defined
%              as zero at the local horizon and increases to pi/2 at
%              zenith.
%
%    rotRad -- Rotation angle of the quad's body frame about the unit vector
%              pointing from the secondary to the primary antenna (i.e.,
%              toward the camera), in rad.  
%
%        rP -- 3x1 position of the quad's primary GNSS antenna, in ECEF
%              coordinates relative to the reference antenna, in meters.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
% OUTPUTS
%
% rCL -------- 3x1 position of the camera in the local ENU frame, whose
%              origin is P.sensorParams.r0G.
%
% RCL -------- 3x3 attitude matrix between the camera and local frames
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Todd Humphreys
%+==============================================================================+  

% Convert az, el, and rotation angles to a body attitude matrix using a 3-2-1
% Euler angle conversion.  3-2-1 is the right sequence because the azimuth
% acts first to change yaw, then elevation acts to change pitch, then roll
% acts.  We subtract pi/2 from azRad because the ENU frame is such that RBL
% should be the identity matrix when azRad = pi/2, elRad = 0, and rotRad = 0.
azRad = x(1); elRad = x(2); rotRad = x(3);
% Rotation matrix that converts a 3x1 coordinate vector expressed in the local
% ENU frame to one expressed in the body frame: vB = RBL*vL
RBL = euler2dcm_321([rotRad;elRad;azRad - pi/2]);

% Form camera-to-local attitude matrix
RCL = P.sensorParams.RCB*RBL;

% Put rP in local coordinates
rP = x(4:6);
RLG = Recef2enu(P.sensorParams.r0G);
rPL = RLG*(rP + P.sensorParams.rRefG - P.sensorParams.r0G);

% Determine location of camera in local coordinates
rCL = rPL + RBL'*(P.sensorParams.rc - P.sensorParams.rA(:,1)); 
