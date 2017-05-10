function [] = PlotEpipoles
% Given two images, their calibrations, correpondences, and the fundamental
% matrix between them plot the epipoles and epipolar lines on the images.
% Inputs:
% img1 and img2: the raw images
% K1, K2: (3x3) calibration matrices for each camera
% k1, k2: (kx1) distortion parameters for each camera
% corrs1, corrs2: (nx2) lists of corresponding points in each image
% F: (3x3) fundamental matrix mapping from First camera's frame to second
% camera's frame

load('my data/cam1.mat');
load('my data/cam2.mat');
load('my data/correspondences3_afterransac.mat');
% F = EstimateFundamentalMatrix(corrs1,corrs2); % for initial estimated F
[ ~, F, ~, ~, ~ ] = GetTransformation( corrs1, corrs2, K1, K2 ); % for reconstucted F
img1 = imread('my data/im3(1)_cam1.jpg');
img2 = imread('my data/im3(1)_cam2.jpg');

% load('test_images/intrinsic.mat');
% K1=K;
% K2=K;
% k1=[0;0];
% k2=[0;0];
% load('test_images/correspondence_ransac0.005.mat');
% % F = EstimateFundamentalMatrix(corrs1,corrs2); % for initial estimated F
% [ ~, F, ~, ~, ~ ] = GetTransformation( corrs1, corrs2, K1, K2 ); % for reconstucted F
% img1 = imread('test_images/l_undistorted.bmp');
% img2 = imread('test_images/r_undistorted.bmp');


img1 = myUndistortImage(img1, K1, k1);
img2 = myUndistortImage(img2, K2, k2);

%% TODO for student: Implement myEpipolarLine and ComputeEpipole
lines1 = myEpipolarLine(F',corrs2);
lines2 = myEpipolarLine(F, corrs1);

[ep1, ep2] = ComputeEpipole(F);

%% Visualization
% Display Image 1, and the distorted Image 2
f0 = figure();
f1 = subplot(1,2,1); imshow(img1);
f2 = subplot(1,2,2); imshow(img2);

% draw epipolar line on left image
set(f0, 'currentaxes', f1), hold on
x = [1, size(img1,2)];
imshow(img1);
ax = gca;

% Plot the lines and points on the left inage
for j = 1 : size(lines1, 1)
    y = -(lines1(j,1) * x + lines1(j,3)) / lines1(j,2);
    plot(x, y);
    ax.ColorOrderIndex = max(1, ax.ColorOrderIndex - 1);
    plot(corrs1(j,1), corrs1(j,2), 'x');
end
% Giant black and white x marks the epipole!
plot(ep1(1), ep1(2), 'wx' , 'markers',18, 'linewidth',2);
plot(ep1(1), ep1(2), 'kx', 'markers',12);
hold off

% Plot the lines and points on the right inage
set(f0, 'currentaxes', f2), hold on
x = 0 : 0.1 : size(img2,2);
imshow(img2);
ax = gca;

for j = 1 : size(lines2, 1)
    y = -(lines2(j,1) * x + lines2(j,3)) / lines2(j,2);
    plot(x, y);
    ax.ColorOrderIndex = max(1, ax.ColorOrderIndex - 1);
    plot(corrs2(j,1), corrs2(j,2), 'x');
end
plot(ep2(1), ep2(2), 'wx' , 'markers',18, 'linewidth',2);
plot(ep2(1), ep2(2), 'kx' , 'markers',12);
hold off

end

