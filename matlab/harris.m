function [ x, y, scores, Ix, Iy ] = harris( image )
% Credit: https://tonypoer.io/2016/10/01/experimenting-with-the-harris-corner-detector-algorithm-in-matlab/
% Work originally by Anthony Poerio
% Modifications made:
%   - remove Gaussian blur on image gradients, not necessary in this case
%   - update thresholds to fit with test images

% HARRIS_CORNERS Extracts points with a high degree of 'cornerness' from
% RGB image matrix of type uint8
%   Input - image = NxMx3 RGB image matrix
%   Output - x = nx1 vector denoting the x location of each of n
%                detected keypoints 
%            y = nx1 vector denoting the y location of each of n 
%                detected keypoints
%            scores = an nx1 vector that contains the value (R) to which a
%            a threshold was applied, for each keypoint
%            Ix = A matrix with the same number of rows and columns as the
%            input image, storing the gradients in the x-direction at each
%            pixel
%            Iy = A matrix with the same nuimber of rwos and columns as the
%            input image, storing the gradients in the y-direction at each
%            pixel

    % compute the gradients, re-use code from HW2P, use window size of 5px
    % convert image to grayscale first
    G = image;

    % convert to double
    G2 = im2double(G);
    
    % create X and Y Sobel filters
    horizontal_filter = [1 0 -1; 2 0 -2; 1 0 -1];
    vertical_filter = [1 2 1; 0 0 0 ; -1 -2 -1];

    % using imfilter to get our gradient in each direction
    filtered_x = conv2(horizontal_filter, G2);
    filtered_y = conv2(vertical_filter, G2);
    
    % store the values in our output variables, for clarity
    Ix = filtered_x;
    Iy = filtered_y;
    
    % Compute the values we need for the matrix...
    % Using a gaussian blur, because I get more positive values after applying
    % it, my values all skew negative for some reason...
    Ix2 = Ix.^2; 
    Iy2 = Iy.^2; 
    Ixy = Ix.*Iy;
    
    % set empirical constant between 0.04-0.06
    k = 0.04;

    num_rows = size(image,1);
    num_cols = size(image,2);

    % create a matrix to hold the Harris values
    H = zeros(num_rows, num_cols);

    % get our matrix M for each pixel
    for y = 6:size(image,1)-6         % avoid edges
        for x = 6:size(image,2)-6     % avoid edges  
            % calculate means (because mean is sum/num pixels)
            % generally, this algorithm calls for just finding a sum,
            % but using the mean makes visualization easier, in my code,
            % and it doesn't change which points are computed to be corners.
            % Ix2 mean
            Ix2_matrix = Ix2(y-2:y+2,x-2:x+2);
            Ix2_mean = sum(Ix2_matrix(:));

            % Iy2 mean
            Iy2_matrix = Iy2(y-2:y+2,x-2:x+2);
            Iy2_mean = sum(Iy2_matrix(:));

            % Ixy mean
            Ixy_matrix = Ixy(y-2:y+2,x-2:x+2);
            Ixy_mean = sum(Ixy_matrix(:));

            % compute R, using te matrix we just created
            Matrix = [Ix2_mean, Ixy_mean; 
                      Ixy_mean, Iy2_mean];
            R1 = det(Matrix) - (k * trace(Matrix)^2);
            
            % store the R values in our Harris Matrix
            H(y,x) = R1;

        end
    end
    
    % set threshold of 'cornerness' 
    avg_r = mean(mean(H));
    threshold = abs(4 * avg_r);

    [row, col] = find(H > threshold);

    scores = [];
    %get all the values
    for index = 1:size(row,1)
        %see what the values are
        r = row(index);
        c = col(index);
        scores = cat(2, scores,H(r,c));
    end

    y = row;
    x = col;

end