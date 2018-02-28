function [Ib] = gaussian_blur(I, wndSize, sigma)
% GAUSSIAN_BLUR Smooth image with symmetric Gaussian filter.
%
%   [Ib] = GAUSSIAN_BLUR(I, wndSize, sigma) produces a filtered image Ib 
%   from I using a square Gaussian kernel with window size wndSize.
%
%   Inputs:
%   -------
%    I        - mxn intensity image.
%    wndSize  - Kernel window size (square, odd number of pixels).
%    sigma    - Standard deviation of Gaussian (pixels, symmetric).
%
%   Outputs:
%   --------
%    Ib  - mxn filtered output image, of same size and class as I.

    % Get distance from center for all points in kernel
    halfway = floor(wndSize/2);
    [x, y] = meshgrid(-halfway:halfway, -halfway:halfway);

    % Create Gaussian filter of kernel size
	filter = exp(-(x .^ 2 + y .^ 2) / (2 * sigma ^ 2));
    filter = filter / sum(sum(filter));

    % Perform convolution on image with Gaussian filter and cast back to
    % class I
    Ib = conv2(double(I), filter, 'same');
    Ib = cast(Ib, 'like', I);
  
end
