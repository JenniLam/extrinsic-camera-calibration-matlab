function [rpy] = rpy_from_dcm(R)
% Credit: Code provided by Uoft ROB501 Course by Jonathan Kelly

% RPY_FROM_DCM Roll, pitch, yaw Euler angles from rotation matrix.
%
%   [rpy] = RPY_FROM_DCM(R) computes roll, pitch and yaw angles from the
%   rotation matrix R.  The pitch angle p is constrained to the range
%   (-pi/2, pi/2].  The returned angles are in radians.
%
%   Inputs:
%   -------
%    R  - 3x3 orthonormal rotation matrix.
%
%   Outputs:
%   --------
%    rpy  - 3x1 vector of roll, pitch, yaw Euler angles.

rpy = zeros(3, 1);

% Roll.
rpy(1) = atan2(R(3, 2), R(3, 3));

% Pitch.
sp = -R(3, 1);
cp = sqrt(R(1, 1).*R(1, 1) + R(2, 1).*R(2, 1));

if abs(cp) > 1e-13  % Or maybe this should be exactly zero?
  rpy(2) = atan2(sp, cp);
else
  rpy(2) = pi/2;
  
  if sp < 0
    rpy(2) = -rpy(2);
  end
end

% Yaw.
rpy(3) = atan2(R(2, 1), R(1, 1));