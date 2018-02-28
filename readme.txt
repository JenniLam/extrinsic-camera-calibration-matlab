Extrinsic Camera Calibration

Overview
------------------------------------------------------------------
This code requires an initial guess of the camera location, the
intrinsic camera matrix, images of a checkerboard calibration
target, the 3D world frame coordinates of the checkerboard cross
junctions and the four corners bounding the calibration target in
the image. It then outputs a corrected pose of the camera.

Note:
- the bounding corners of the calibration target do not have to be 
  particularly exact, but are necessary in cases where the
  checkerboard is heavily warped in image
- the images are expected to have little radial or tangential
  distortion

Automatic Cross Junction Detection (cross points in checkerboard):
1. Use four corner bounding box to compute transform from
   warped to flat checkerboard.
2. Given homography, bilinearly interpolate warped pixels and
   assign to locations in flat checkerboard image.
3. Blur image to reduce noise, and run harris detection to find
   corners.
4. Collapse clusters of detected features into centroid of clusters.
   (Harris detector typically finds a cluster of 40 'corners' 
   centered around each cross junction)
5. Run saddle point detection on small patches around detected
   feature cluster centroid. This provides a more exact location for
   the cross junction.
6. Filter out L and T junctions from X junctions.
7. Sort and order cross junctions in row-major form (necessary for
   extrinsic calibration step).
8. Transform back to warped image.

Extrinsic Calibration:
1. Given initial guess of parameters (roll, pitch, yaw, x, y, z),
   estimated cross junctions in images, and actual 3D location of
   cross junctions on the checkerboard, compute Jacobian and
   error residuals.
2. Use Gauss-Newton to improve parameters, and converge on best
   guess for camera pose.

Credit
------------------------------------------------------------------
Uses code originally by Anthony Poerio for Harris corner detection
https://tonypoer.io/2016/10/01/experimenting-with-the-harris-
corner-detector-algorithm-in-matlab/
See harris.m

Uses code provided by University of Toronto ROB501 Course by 
Jonathan Kelly for transformations between rotation matrices and
roll, pitch, yaw
See dcm_from_rpy.m, dcm_jacob_rpy.m, rpy_from_dcm.m

Implements cross junction and saddle point detection originally
described by:
L. Lucchese and S. K. Mitra, "Using Saddle Points for Subpixel Feature
Detection in Camera Calibration Targets," in Proc. Asia-Pacific Conf.
Circuits and Systems (APCCAS'02), vol. 2, (Singapore), pp. 191-195,
Dec. 2002.
See saddle_point.m

Contents
-------------------------------------------------------------------
images/                     images tested on
matlab/                     MATLAB code
test/                       test files, including ground truth data

To Run
-------------------------------------------------------------------
Run run.m
(Tested with MATLAB 2015)