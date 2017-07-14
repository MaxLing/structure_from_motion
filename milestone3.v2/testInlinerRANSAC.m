img1 = 1;
img2 = 2;

%% Data Management
K = [568.996140852 0 643.21055941;
     0 568.988362396 477.982801038;
     0 0 1];
img_num = 6;

% Load images
for i = 1 : img_num
    str = sprintf('Milestone3_data/SfMProjectData_1/image%07d.bmp', i);
    img{i} = imread(str);
end

% load matching
MatchX = []; % n*6
MatchY = []; % n*6
Match  = []; % n*6 n*15
for i = 1 : img_num-1;
    str = sprintf('Milestone3_data/SfMProjectData_1/matching%d.txt', i);
    [matchX, matchY, match] = LoadMatching(str, i, img_num);
    MatchX = [MatchX;matchX];
    MatchY = [MatchY;matchY];
    Match = [Match;match];
end

% initialize 3d points and visibility matrix
X3D = zeros(size(MatchX,1), 3);
V = zeros(size(MatchX,1), img_num);

idx1 = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-img1))+img2-img1)==1);
x1 = [MatchX(idx1,img1) MatchY(idx1,img1)];
x2 = [MatchX(idx1,img2) MatchY(idx1,img2)];
figure;
showMatchedFeatures(img{img1},img{img2},x1,x2, 'method', 'montage');

%%  Matching RANSAC
        idx = find(Match(:,img1)==1 & Match(:,img2)==1); % feature matching in two images
        x1 = [MatchX(idx,img1) MatchY(idx,img1)];
        x2 = [MatchX(idx,img2) MatchY(idx,img2)];
        [x1_new,x2_new,inliners_id] = GetInlinersRANSAC(x1,x2);
        Match(setxor(idx,idx(inliners_id)),img1) = 0; % remove outliers

figure;
showMatchedFeatures(img{img1},img{img2},x1_new,x2_new, 'method', 'montage');