%% Load images and data 

% 5 given test images
I1 = imread('../images/target_01.png');
I2 = imread('../images/target_02.png');
I3 = imread('../images/target_03.png');
I4 = imread('../images/target_04.png');
I5 = imread('../images/target_05.png');

% bounds for the checkerboard in each image
load('../test/bounds.mat');

% ground truth cross junction points in the image
% used to evaluate results
load('../test/cj_pts_image_01.mat');
cj_pts1 = cj_pts;
load('../test/cj_pts_image_02.mat');
cj_pts2 = cj_pts;
load('../test/cj_pts_image_03.mat');
cj_pts3 = cj_pts;

% actual location of the cross junctions of the checkerboard in the world
% frame, required for pose estimation of camera
load('../test/Wpts.mat');

% load the initial guess of the camera location 'Eg', and the ground truth
% camera position 'Ef' for evaluating camera calibration
load('../test/camera_pose_01.mat');
Eg1 = [camera.guess.C camera.guess.t];
Eg1 = [Eg1; 0 0 0 1];
Ef1 = [camera.best.C  camera.best.t];
Ef1 = [Ef1; 0 0 0 1];
load('../test/camera_pose_02.mat');
Eg2 = [camera.guess.C camera.guess.t];
Eg2 = [Eg2; 0 0 0 1];
Ef2 = [camera.best.C  camera.best.t];
Ef2 = [Ef2; 0 0 0 1];
load('../test/camera_pose_03.mat');
Eg3 = [camera.guess.C camera.guess.t];
Eg3 = [Eg3; 0 0 0 1];
Ef3 = [camera.best.C  camera.best.t];
Ef3 = [Ef3; 0 0 0 1];
load('../test/camera_pose_04.mat');
Eg4 = [camera.guess.C camera.guess.t];
Eg4 = [Eg4; 0 0 0 1];
Ef4 = [camera.best.C  camera.best.t];
Ef4 = [Ef4; 0 0 0 1];
load('../test/camera_pose_05.mat');
Eg5 = [camera.guess.C camera.guess.t];
Eg5 = [Eg5; 0 0 0 1];
Ef5 = [camera.best.C  camera.best.t];
Ef5 = [Ef5; 0 0 0 1];

%% Find cross junctions on checkerboard

% find cross junctions
pts1 = cross_junctions(I1, bpoly1, Wpts);
pts2 = cross_junctions(I2, bpoly2, Wpts);
pts3 = cross_junctions(I3, bpoly3, Wpts);
pts4 = cross_junctions(I4, bpoly4, Wpts);
pts5 = cross_junctions(I5, bpoly5, Wpts);

% display found cross junctions
figure(1);
hold off;
imshow(I1);
hold on;
plot(pts1(1,:), pts1(2,:), 'x');
k = 1:size(pts1, 2);
text(pts1(1,:), pts1(2,:), num2str(k'), 'Color', 'r');

figure(2);
hold off;
imshow(I2);
hold on;
plot(pts2(1,:), pts2(2,:), 'x');
k = 1:size(pts2, 2);
text(pts2(1,:), pts2(2,:), num2str(k'), 'Color', 'r');

figure(3);
hold off;
imshow(I3);
hold on;
plot(pts3(1,:), pts3(2,:), 'x');
k = 1:size(pts3, 2);
text(pts3(1,:), pts3(2,:), num2str(k'), 'Color', 'r');

figure(4);
hold off;
imshow(I4);
hold on;
plot(pts4(1,:), pts4(2,:), 'x');
k = 1:size(pts4, 2);
text(pts4(1,:), pts4(2,:), num2str(k'), 'Color', 'r');

figure(5);
hold off;
imshow(I5);
hold on;
plot(pts5(1,:), pts5(2,:), 'x');
k = 1:size(pts5, 2);
text(pts5(1,:), pts5(2,:), num2str(k'), 'Color', 'r');

% Compute errors on found cross junctions
fprintf('Image 1 cross junction errors: %f\n', sum(sum(abs(cj_pts1 - pts1)))/size(cj_pts1,2));
fprintf('Image 2 cross junction errors: %f\n', sum(sum(abs(cj_pts2 - pts2)))/size(cj_pts2,2));
fprintf('Image 3 cross junction errors: %f\n\n', sum(sum(abs(cj_pts3 - pts3)))/size(cj_pts3,2));

%% Run pose estimation

% run pose estimation
E1 = pose_estimate_nlopt(Eg1, pts1, Wpts);
E2 = pose_estimate_nlopt(Eg2, pts2, Wpts);
E3 = pose_estimate_nlopt(Eg3, pts3, Wpts);
E4 = pose_estimate_nlopt(Eg4, pts4, Wpts);
E5 = pose_estimate_nlopt(Eg5, pts5, Wpts);

% compute and display errors on pose estimation
aa1 = rotm2axang(E1(1:3,1:3)'*Ef1(1:3,1:3));
aa2 = rotm2axang(E2(1:3,1:3)'*Ef2(1:3,1:3));
aa3 = rotm2axang(E3(1:3,1:3)'*Ef3(1:3,1:3));
aa4 = rotm2axang(E4(1:3,1:3)'*Ef4(1:3,1:3));
aa5 = rotm2axang(E5(1:3,1:3)'*Ef5(1:3,1:3));

fprintf('Image 1 translational error: %f, rotational error %f\n', norm(E1(1:3,4)-Ef1(1:3,4)), norm(aa1(1:3)*aa1(4)/norm(aa1(1:3))));
fprintf('Image 2 translational error: %f, rotational error %f\n', norm(E2(1:3,4)-Ef2(1:3,4)), norm(aa2(1:3)*aa2(4)/norm(aa2(1:3))));
fprintf('Image 3 translational error: %f, rotational error %f\n', norm(E3(1:3,4)-Ef3(1:3,4)), norm(aa3(1:3)*aa3(4)/norm(aa3(1:3))));
fprintf('Image 4 translational error: %f, rotational error %f\n', norm(E4(1:3,4)-Ef4(1:3,4)), norm(aa4(1:3)*aa4(4)/norm(aa4(1:3))));
fprintf('Image 5 translational error: %f, rotational error %f\n', norm(E5(1:3,4)-Ef5(1:3,4)), norm(aa5(1:3)*aa5(4)/norm(aa5(1:3))));