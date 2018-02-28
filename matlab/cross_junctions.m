function [Ipts] = cross_junctions(I, boundPoly, Wpts)
% CROSS_JUNCTIONS Find cross-junctions in image with subpixel accuracy.
%
%   [Ipts] = CROSS_JUNCTION(I, boundPoly, Wpts) locates a series of cross- 
%   junction points on a planar calibration target, where the target is
%   bounded in the image by the specified 4-sided polygon. The number of
%   cross-junctions identified should be equal to the number of world points.
%
%   Note also that the world and image points must be in *correspondence*,
%   that is, the first world point should map to the first image point, etc.
%
%   Inputs:
%   -------
%    I          - Image (grayscale, double or integer class).
%    boundPoly  - 2x4 array defining bounding polygon for target (clockwise).
%    Wpts       - 3xn array of world points (in 3D, on calibration target).
%
%   Outputs:
%   --------
%    Ipts  - 2xn array of cross-junctions (x, y), relative to the upper left
%            corner of I.

% STEP 1. Pre-process image
%   - gaussian blur image in preparation for corner detection and 
%     saddle point location
%   - transform checkerboard from warped -> flat

    % From dimensions of board, get corners of transformed, flat checkerboard
    board  = [1 572 572 1;
              1 1    445 445];

    % Sort given corners of board in image
    [~, idxX] = sort(boundPoly(1,:));
    [~, idxY] = min(boundPoly(2,idxX(1:2)));
    
    % If corners of board are not in order from (top,left) -> (bot,left)
    % (like dimensions of board), shift them into this order
    % This is necessary to match correspondences of homography
    if     idxX(idxY) == 1; numShift = 0;
    elseif idxX(idxY) == 2; numShift = 3;
    elseif idxX(idxY) == 3; numShift = 2;
    else                    numShift = 1;
    end
    boundPoly = circshift(boundPoly, numShift, 2);
    
    % Compute homography from flat -> warped
    H = dlt_homography(board, boundPoly);
    
    % Generate all pixel locations in flat checkerboard
    [hor, ver] = meshgrid(1:572,1:445);
    BoardPts = [hor(:)'; ver(:)'; ones(1,length(hor(:)))];
    
    % Transform all pixels from flat to warped checkerboard
    Bpts = H * BoardPts;
    Bpts(1,:) = Bpts(1,:) ./ Bpts(3,:);
    Bpts(2,:) = Bpts(2,:) ./ Bpts(3,:);
    
    % Get new image by bilinearly interpolating pixel intensities from
    % warped pixel locations and assigning to flat checkerboard image
    Im    = zeros(445, 572);
    Im(:) = bilinear_interp(I, Bpts);
    [height, width] = size(Im);
    
% STEP 2. Get Preliminary Saddle Point Locations
%   - Run Harris corner detector to gather all corner locations
%   - From output of corner detection, narrow down preliminary clusters
%     of corner locations to single points
    
    % Blur image
    Im = gaussian_blur(Im, 3, 1);

    % Extract x,y locations of corners through Harris corner detector
    % See original source code: https://tonypoer.io/2016/10/01/experimenting-with-the-harris-corner-detector-algorithm-in-matlab/
    [x, y, ~, ~, ~] = harris(Im);
    pts = [x'; y'];

    % Set maxdist, the maximum distance 2 pixels can be apart, and still be
    % considered to be in the same 'cluster'
    single_pts = [];
    maxdist    = 10;
    
    % Gather 'cluster' of points (all points within 10 pixels of each
    % other)
    % Take average of all point locations in cluster, use this as the
    % preliminary X-Junction location
    while size(pts,2) > 0
        
        % Get all points within 10 pixels 
        group   = sqrt((pts(1,:)-pts(1,1)).^2 + (pts(2,:)-pts(2,1)).^2) < maxdist;
        found   = pts(:, group);
        
        % Get all points within 10 pixels of all the points in the cluster,
        % and all points within 10 pixels of these points, and so
        % on...
        i = 1;
        while i <= size(found,2)
            nextpts = sqrt((pts(1,:)-found(1,i)).^2 + (pts(2,:)-found(2,i)).^2) < maxdist;
            nextpts = ~group & nextpts;
            group   = group | nextpts;
            found   = [found pts(:, nextpts)];
            i = i+1;
        end
        
        % Take the average location of the points in the cluster, and add
        % it as the location which represents the cluster
        grouped_pts = pts(:, group);
        avg_pt      = [round(mean(grouped_pts(1,:))); round(mean(grouped_pts(2,:)))];
        single_pts  = [single_pts avg_pt];
        
        % Remove all points in the cluster from the set of points (avoids
        % extra computations and duplicate points)
        pts(:, group) = [];
    end
    
% STEP 3. Run Saddle Point Locator
%   - For the single points detected, create patches from surrounding
%     pixels, prune out L junctions in checkerboard, and run saddle point
%     locator on patches
    
    % Set the half width of the patch and the XJunctionScale, used to
    % differentiate between X and L junctions
    halfPatchWidth = 5;
    XJunctionScale = 1.3;
    
    % Iterate through all single points generated in previous steps, and
    % get specific saddle point locations
    for k = 1:3
        
        unsorted_pts = zeros(2, size(Wpts, 2));
        fidx = 1;
        
        for i = 1:size(single_pts,2)

            % Get corners of patch, based on chosen halfPatchWidth
            left  = min(max(single_pts(1,i)-halfPatchWidth, 1), width-2*halfPatchWidth);
            right = left + 2*halfPatchWidth;
            top   = min(max(single_pts(2,i)-halfPatchWidth, 1), height-2*halfPatchWidth);
            bot   = top + 2*halfPatchWidth;

            % Compute total intensity of N, S, W, E borders of patch (use 3
            % pixel thick border)
            % Sort these total intensities to compare X vs. L junction
            borders = [ sum(sum(Im(top:(top+3), left:right)));
                        sum(sum(Im((bot-3):bot, left:right)));
                        sum(sum(Im(top:bot,     left:(left+3))));
                        sum(sum(Im(top:bot,     (right-3):right)))];
            borders = sort(borders);

            % If patch is an X junction, run saddle point detection
            if (borders(1)+borders(2))*XJunctionScale > (borders(3)+borders(4))
                spt = saddle_point(Im(top:bot, left:right));

                % If saddle point is found within patch, add to array of
                % X-junction points
                if spt(1) <= 2*halfPatchWidth && spt(1) >=1 && spt(2) <= 2*halfPatchWidth && spt(2) >= 1
                    unsorted_pts(1,fidx) = left + spt(1);
                    unsorted_pts(2,fidx) = top  + spt(2);
                    fidx = fidx + 1;
                end
            end  
        end
        
        % run saddle point detection again and threshold higher on X
        % junctions when not enough points are detected
        if unsorted_pts(1,size(Wpts, 2)) ~= 0 && unsorted_pts(2,size(Wpts, 2)) ~= 0
            break;
        else
           XJunctionScale = XJunctionScale + 1.3;
        end
    end

% STEP 4. Post-process Detected Saddle Points
%   - Sort points in row-major order
%   - Transform back to original locations using previously computed
%     homography
    
    % Sort points based on Y location: first 8 points will be first row,
    % second 8 points will be second row, etc.
    sortedY = sortrows(unsorted_pts', 2);
    
    % Sort points in each row based on X location
    SIpts = zeros(3,48);
    SIpts(1:2,1:8)   = sortrows(sortedY(1:8,:))';
    SIpts(1:2,9:16)  = sortrows(sortedY(9:16,:))';
    SIpts(1:2,17:24) = sortrows(sortedY(17:24,:))';
    SIpts(1:2,25:32) = sortrows(sortedY(25:32,:))';
    SIpts(1:2,33:40) = sortrows(sortedY(33:40,:))';
    SIpts(1:2,41:48) = sortrows(sortedY(41:48,:))';
    SIpts(3,:)       = 1;
    
    % Transform points back to warped checkerboard
    SIpts = H*SIpts;
    Ipts  = [SIpts(1,:)./SIpts(3,:); SIpts(2,:)./SIpts(3,:)];
  
end