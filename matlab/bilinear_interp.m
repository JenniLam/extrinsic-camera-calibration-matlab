function [b] = bilinear_interp(I, pt)
% bilinear_interp Performs bilinear interpolation for a given image point.
%
%   Given the (x y) location of a point in an input image, use the 
%   surrounding 4 pixels to output the bilinearly interpolated intensity.
%
%   Note that images are (usually) integer-valued functions (in 2D), therefore
%   the intensity value you return should be an integer.
%
%  Inputs:
%  -------
%   I   - Input image (monochrome, one channel - n rows x m columns).
%   pt  - Point in input image (x, y), with subpixel precision.
%
%  Outputs
%  -------
%   b  - Interpolated brightness or intensity value (whole number >= 0).

    % Get dimensions of image
    [height, width]  = size(I);
    
    % Get neighbouring pixels
    top   = max(floor(pt(2,:)),1); % top
    bot   = min(top+1,height);   % bottom
    left  = max(floor(pt(1,:)),1); % left
    right = min(left+1,width);   % right
    
    % Perform bilinear interpolation
    idxTR = sub2ind(size(I), top, right);
    idxTL = sub2ind(size(I), top, left);
    idxBR = sub2ind(size(I), bot, right);
    idxBL = sub2ind(size(I), bot, left);
    
    top_b = (pt(1,:)-left).*double(I(idxTR)) + (right-pt(1,:)).*double(I(idxTL));
    bot_b = (pt(1,:)-left).*double(I(idxBR)) + (right-pt(1,:)).*double(I(idxBL));
    b = ceil((pt(2,:)-top).*bot_b + (bot-pt(2,:)).*top_b);

end
