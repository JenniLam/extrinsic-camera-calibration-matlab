function [E] = pose_estimate_nlopt(Eg, Ipts, Wpts)
%  POSE_ESTIMATE_NLOPT Estimate camera pose from 2D-3D correspondences via NLS.
%
%   [E] = POSE_ESTIMATE_NLOPT(Eg, Ipts, Wpts) performs a nonlinear least squares 
%   optimization procedure to determine the best estimate of the camera pose in 
%   the calibration target frame, given 2D-3D point correspondences.
%
%   Inputs:
%   -------
%    Eg    - 4x4 homogenous pose matrix, initial guess for camera pose.
%    Ipts  - 2xn array of cross-junction points (with subpixel accuracy).
%    Wpts  - 3xn array of world points (one-to-one correspondence with image).
%
%   Outputs:
%   --------
%    E  - 4x4 homogenous pose matrix, estimate of camera pose in target frame.

    % Get initial guess
    E = Eg;
    
    % Number of iterations of non-linear least squares
    max_iter = 10;
    
    % Number of points
    num_pts = size(Ipts, 2);
    
    % Camera intrinsic matrix
    K = [564.9 0 337.3; 0 564.3 226.5; 0 0 1];
    
    % Parameters to adjust
    p = zeros(6,1);
    R = E(1:3, 1:3);
    RPY = rpy_from_dcm(R);
    p(4:6) = RPY;
    p(1:3) = E(1:3,4);
    
    % Set up Jacobian
    J = zeros(2*num_pts, 6);
    
    % Set up residuals
    error = zeros(2*num_pts, 1);
    
    iter = 1;
    while iter <= max_iter
        
        % Get Homogeneous transformation from parameters
        E = get_homogeneous_transform(p);
        
        % Calculate Jacobian and residuals for all cross junctions
        for i = 1:size(Ipts, 2)
            J((2*i-1):(2*i),:)   = find_jacobian(K, E, Wpts(:,i));
            error((2*i-1):(2*i)) = transform_point(K, E, Wpts(:,i)) - Ipts(:,i);
        end
        
        % Update parameters
        deltaP = pinv(J)*error;
        p = p - deltaP;
        iter = iter+1;
    end

    % Get Homogeneous transformation from parameters
    function H = get_homogeneous_transform(p)
        Rot = dcm_from_rpy([p(4), p(5), p(6)]);
        H = zeros(4,4);
        H(1:3,1:3) = Rot;
        H(1:3,4)   = [p(1); p(2); p(3)];
        H(4,4)     = 1;
    end

    % Get point in image plane given camera intrinsic matrix, homogeneous 
    % transform,and 3D point
    function tpt = transform_point(K, H, world_pt)
        trans_pt = world_pt - H(1:3,4);
        trans_pt = K*H(1:3,1:3)'*trans_pt;
        tpt      = trans_pt(1:2)/trans_pt(3);
    end

end

 
