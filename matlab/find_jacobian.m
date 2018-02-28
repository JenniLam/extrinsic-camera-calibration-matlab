function [J] = find_jacobian(K, Ec, Wpt)
%  FIND_JACOBIAN Determine Jacobian for NLS camera pose optimization.
%
%   [J] = FIND_JACOBIAN(K, Ec, Wpt) computes the Jacobian of image plane point
%   with respect to the current camera pose estimate and given a single world
%   point.
%
%   Inputs:
%   -------
%    K    - 3x3 camera intrinsic calibration matrix.
%    Ec   - 4x4 homogenous pose matrix, current guess for camera pose.
%    Wpt  - 3x1 world point on calibration target (one of n).
%
%   Outputs:
%   --------
%    J  - 2x6 Jacobian matrix (columns are tx, ty, tz, r, p, q).

    % homogeneous transform and real world point
    Hwc = Ec;
    P = Wpt;

    % rotation matrix
    R  = Hwc(1:3, 1:3);
    
    % difference between world point and origin of camera
    dx = P - Hwc(1:3, 4);

    % separate out last row, since this is the normalization term for the
    % points in the image plane, it ends up as the denominator for x and y
    % in the image plane
    l = K*R'*dx;
    g = l(3);

    dldP = K*R.';

    % translational derivatives
    dldH(1:3, 1:3) = -dldP;

    % get derivatives of rotation matrix with respect to roll, pitch, and
    % yaw and compute final rotational derivatives
    [dRdr, dRdp, dRdq] = dcm_jacob_rpy(R);
    dldH(1:3, 4) = K*dRdr'*dx;
    dldH(1:3, 5) = K*dRdp'*dx;
    dldH(1:3, 6) = K*dRdq'*dx;

    dgdH = dldH(3, :);

    % use chain rule gives final jacobian
    dpdH = (g*dldH - l*dgdH)/(g^2);
    J = dpdH(1:2, :);

end
