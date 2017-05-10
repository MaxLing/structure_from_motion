load('my data/cam1.mat');
load('my data/cam2.mat');

img1 = imread('my data/im3(1)_cam1.jpg');
img2 = imread('my data/im3(1)_cam2.jpg');

% load('test_images/intrinsic.mat');
% K1=K;
% K2=K;
% k1=[0;0];
% k2=[0;0];
% img1 = imread('test_images/l_undistorted.bmp');
% img2 = imread('test_images/r_undistorted.bmp');


img1 = myUndistortImage(img1, K1, k1);
img2 = myUndistortImage(img2, K2, k2);
img1 = rgb2gray(img1);
img2 = rgb2gray(img2);

points1 = detectSURFFeatures(img1);
points2 = detectSURFFeatures(img2);

[features1,valid_points1] = extractFeatures(img1,points1);
[features2,valid_points2] = extractFeatures(img2,points2);

indexPairs = matchFeatures(features1,features2);

corr1 = valid_points1(indexPairs(:,1),:);
corr2 = valid_points2(indexPairs(:,2),:);

% ransac
[corr1_new,corr2_new] = ransac(corr1,corr2,K1,K2,8,0.99,0.005,0.5);

figure;
showMatchedFeatures(img1,img2,corr1,corr2, 'method', 'montage');
figure;
showMatchedFeatures(img1,img2,corr1_new,corr2_new, 'method', 'montage');