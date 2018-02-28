function [pt] = saddle_point(I)
% SADDLE_POINT Locate saddle point in an image patch.
%
%   [pt] = SADDLE_POINT(I) finds the subpixel center of a cross-junction in the 
%   image patch I, by blurring the patch, fitting a hyperbolic paraboloid to 
%   it, and then finding the critical point of that paraboloid.
%
%   Note that the location of 'p' is relative to (0.5, 0.5) at the upper left 
%   corner of the patch, i.e., the pixels are treated as covering an area of
%   one unit square.
%
%   Inputs:
%   -------
%    I  - mxn image patch (grayscale, double or integer class).
%
%   Outputs:
%   --------
%    pt  - 2x1 subpixel location of saddle point in I (x, y coords).
%
% References:
%
%   L. Lucchese and S. K. Mitra, "Using Saddle Points for Subpixel Feature
%   Detection in Camera Calibration Targets," in Proc. Asia-Pacific Conf.
%   Circuits and Systems (APCCAS'02), vol. 2, (Singapore), pp. 191-195,
%   Dec. 2002.

    % Image has been pre-blurred
    [height, width] = size(I);
    
    % Get all points in patch
    [hor, ver] = meshgrid(1:width,1:height);
    pts = [ver(:) hor(:)];
    
    % Generate rows of [1 y x y^2 xy x^2] for each point in the patch
    X = zeros(height*width, 6);
    X(:,1)   = 1;
    X(:,2:3) = pts;
    X(:,4)   = X(:,2).^2;
    X(:,5)   = X(:,2).*X(:,3);
    X(:,6)   = X(:,3).^2;
    
    % Generate a Z vector of intensities at all points in the patch
    Z = zeros(height*width, 1);
    for i = 1:height*width
        Z(i) = double(I(X(i,2), X(i,3)));
    end
    
    % Solve for a hyperboloid which fits these intensities
    param = pinv(X)*Z;
    
    % Find saddle point given the fit hyperboloid
    pt = -inv([2*param(6) param(5); param(5) 2*param(4)])*[param(3); param(2)];
    
end