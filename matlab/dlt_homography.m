function [H, A] = dlt_homography(I1pts, I2pts)
% dlt_homography Perspective Homography between two images.
%
%   Given 4 points from 2 separate images, compute the perspective homograpy
%   (warp) between these points using the DLT algorithm.
%
%   Inputs:
%   -------
%    I1pts  - 2x4 array of points from Image 1 (each column is x, y).
%    I2pts  - 2x4 array of points from Image 2 (1-to-1 correspondence).
%
%   Outputs:
%   --------
%    H  - 3x3 perspective homography (matrix map) between image coordinates.
%    A  - 8x9 DLT matrix used to determine homography.

    % matrix map and DLT matrix
    A = zeros(8,9);
    H = zeros(3,3);

    % Compute similarity Transform T1 for I1pts
    T1 = compute_similarity_transform(I1pts);

    % Compute similarity Transform T2 for I2pts
    T2 = compute_similarity_transform(I2pts);
    
    % Get DLT matrix A based on transformed I1pts and I2pts
    A(1:2,:) = generate_matrix_from_correspondences(I1pts(:,1), I2pts(:,1), T1, T2);
    A(3:4,:) = generate_matrix_from_correspondences(I1pts(:,2), I2pts(:,2), T1, T2);
    A(5:6,:) = generate_matrix_from_correspondences(I1pts(:,3), I2pts(:,3), T1, T2);
    A(7:8,:) = generate_matrix_from_correspondences(I1pts(:,4), I2pts(:,4), T1, T2);
    
    % Get null space of A (solution space for H vector)
    h = null(A);
    if ~any(h)
       error('Cannot find solution for DLT homography!'); 
    end
    H = [h(1) h(2) h(3); h(4) h(5) h(6); h(7) h(8) h(9)];
    
    % Transform H back to given coordinates
    H = inv(T2) * H * T1;
  
    % Given set of points, compute similarity transform such that their origin
    % is the centroid of all points, and their average distance from origin is
    % sqrt(2)
    function T = compute_similarity_transform(pts)

        % Get location of centroid
        avg = mean(pts,1);
        avgX = avg(1);
        avgY = avg(2);

        % Get average distance to centroid of all points
        avgDist = 0;
        for i = 1:4
            avgDist = avgDist + sqrt((pts(1,i)-avgX)^2 + (pts(2,i)-avgY)^2);
        end

        % Using average distance to centroid, calculate required scale factor
        scale = sqrt(2) / (avgDist/4);

        % Using location of centroid and required scale factor, calculate
        % similarity transform
        T = eye(3,3);
        T(1:2,1:2) =  scale * eye(2,2);
        T(1:2,3)   = -scale * [avgX; avgY];
    end

    % Generate rows for DLT matrix given points and similarity transforms
    function a = generate_matrix_from_correspondences(p1, p2, T1, T2)
      % Transform points given similarity transforms
      p1 = T1*[p1; 1];
      p2 = T2*[p2; 1];
      p1 = p1 / p1(3);
      p2 = p2 / p2(3);

      % Generate rows with transformed points
      a = [ -p1(1) -p1(2) -1 0      0       0 p2(1)*p1(1) p2(1)*p1(2) p2(1);
            0      0      0  -p1(1) -p1(2) -1 p2(2)*p1(1) p2(2)*p1(2) p2(2)];
    end
  
end


