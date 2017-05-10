% Run this script for each pair of images.
img1 = imread('image/im1(1)_cam1.jpg');
img2 = imread('image/im1(1)_cam2.jpg');

K1 = [2.768201922155784e+03,-3.283946709730422,1.638682157182430e+03;0,2.767657671862577e+03,1.244945553677645e+03;0,0,1]; 
% Should be 3 x 3, calibration matrix for first camera
K2 = [3.456597718123323e+03, 4.948208007856249,2.009097474856611e+03;0,3.465840808785538e+03,1.535384201592209e+03;0,0,1]; 
% Should be 3 x 3, calibration matrix for second camera

k1 = [0.083541099963479;-0.453999314312728]; % Should be 2 x 1, radial distortion for first camera
k2 = [0.036361182076967; 0.114666012903610]; % Should be 2 x 1, radial distortion for second camera

%% Undistort images
img1_undist = myUndistortImage(img1, K1, k1);
img2_undist = myUndistortImage(img2, K2, k2);

[corr1, corr2] = clickCorrespondences(img1_undist, img2_undist);

%% Visualize points
showMatchedFeatures(img1_undist, img2_undist, corr1, corr2, 'method', 'montage');

% You'll want to modify this to save a separate .mat for each pair of images.
save corresondences1.mat corr1 corr2
