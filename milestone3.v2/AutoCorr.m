K = [2768.201922155784	-3.283946709730422	1638.68215718243;
0.0	2767.657671862577	1244.945553677645;
0.0	0.0	1.0]; % Should be 3 x 3, calibration matrix for first camera

k = [0.083541099963479;
-0.453999314312728]; % Should be 2 x 1, radial distortion for first camera
img_num = 6;

MatchX = []; % n*6
MatchY = []; % n*6
Match  = []; % n*6 n*15
for i = 1:img_num-1
    for j = i+1:img_num
        str1 = sprintf('My Photo/image%d.jpg', i);
        strj = sprintf('My Photo/image%d.jpg', j);
img1 = imread(str1);
img2 = imread(str2);
img1_undist = myUndistortImage(img1, K, k);
img2_undist = myUndistortImage(img2, K, k);
img1_gray = rgb2gray(img1_undist);
img2_gray = rgb2gray(img2_undist);

points1 = detectSURFFeatures(img1_gray);
points2 = detectSURFFeatures(img2_gray);
[features1,valid_points1] = extractFeatures(img1_gray,points1);
[features2,valid_points2] = extractFeatures(img2_gray,points2);
indexPairs = matchFeatures(features1,features2);
corr1 = valid_points1(indexPairs(:,1),:);
corr2 = valid_points2(indexPairs(:,2),:);

