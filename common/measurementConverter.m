function [C] = measurementConverter(M,P)
% measurementConverter : Nonlinear conversion from quad-produced measurements
%                        to a full camera pose measurement expressed in the
%                        local ENU frame whose origin is P.sensorParams.r0G
%
% INPUTS
%
% M ---------- Structure of measurements produced by the quad, with the
%              following elements:
%
%         azRad = Azimuth angle for the unit vector pointing from the
%                 secondary to the primary antenna (i.e., toward the camera),
%                 relative to the local ENU frame, in rad.  Azimuth is defined
%                 as zero toward local north and increases clockwise.
%
%         elRad = Elevation angle for the unit vector pointing from the
%                 secondary to the primary antenna (i.e., toward the camera),
%                 relative to the local ENU frame, in rad.  Elevation is
%                 defined as zero at the local horizon and increases to pi/2
%                 at zenith.
%
%            rP = 3x1 position of the quad's primary GNSS antenna, in ECEF
%                 coordinates relative to the reference antenna, in meters.
%
%    azSigmaRad = Standard deviation of error in azRad, in radians.
%
%    elSigmaRad = Standard deviation of error in elRad, in radians.
%
%            RP = 3x3 error covariance matrix for rP.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
% OUTPUTS
%
% C ---------- Structure of camera pose measurement based on M, with the
%              following elements: 
%
%          rCL = 3x1 position of the camera in the local ENU frame, whose
%                origin is P.sensorParams.r0G.
%
%          RCL = 3x3 attitude matrix between the camera and local frames: vC =
%                RCL*vL.
%
%            R = 6x6 error covariance matrix for the measurements.  Let
%                RCLtrue be the true local-to-camera attitude matrix, and let
%                ea be a 3x1 vector of error Euler angles such that RCLtrue =
%                dRCL(ea)*RCL, where dRCL(ea) is the DCM formed from the error
%                Euler angle vector ea.  Also let er = rCL - rCLtrue, where
%                rCLtrue is the 3x1 true location of the camera center.
%                Finally, let e = [er; ea].  Then R is defined to be equal to
%                E[e*e'].
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Todd Humphreys
%+==============================================================================+  

%----- Setup
nx = 6;
rotRad = 0;  rotSigmaRad = 5*pi/180;
alphaUKF = 1e-3; betaUKF = 2; kappaUKF = 0; 
lambda_p = alphaUKF^2*(kappaUKF + nx) - nx;
c_p = sqrt(nx+lambda_p);
Wm0_p = lambda_p/(nx + lambda_p);
Wmi_p = 1/(2*(nx + lambda_p));
Wc0_p = Wm0_p + 1 - alphaUKF^2 + betaUKF; Wci_p = Wmi_p;

%----- Arrange covariances
xBar = [M.azRad;M.elRad;rotRad;M.rP];
Paer = diag([M.azSigmaRad^2,M.elSigmaRad^2,rotSigmaRad^2]);
Px = blkdiag(Paer,M.RP);
Sx = chol(Px)';

%----- Assemble sigma points and push these through the nonlinear function
sp0 = xBar;
yMat = zeros(nx,2*nx);
[C.rCL,C.RCL] = g_fun(sp0,P);
y0 = [C.rCL;zeros(3,1)];
for ii=1:2*nx
  jj = ii; pm = 1;
  if(ii > nx) jj = ii - nx; pm = -1; end
  spii = sp0 + pm*c_p*Sx(:,jj);
  [rCL,RCL] = g_fun(spii,P);
  yMat(:,ii) = [rCL;dcm2euler(RCL*C.RCL')];
end
% Recombine sigma points
yBar = sum([Wm0_p*y0, Wmi_p*yMat],2);
Py = Wc0_p*(y0 - yBar)*(y0 - yBar)';
for ii=1:2*nx
  Py = Py + Wci_p*(yMat(:,ii) - yBar)*(yMat(:,ii) - yBar)';
end
C.R = Py;