%% Data Management
K = [568.996140852 0 643.21055941;
     0 568.988362396 477.982801038;
     0 0 1];
img_num = 6;

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

%%  Matching RANSAC
for k1 = 1:img_num-1 
    for k2 = k1+1:img_num      % all possible pairs of images
%         idx = find(Match(:,k1)==1 & Match(:,k2)==1); % feature matching in two images
        idx = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-k1))+k2-k1)==1);
        if length(idx) < 8 % not enough correspondence
            continue;
        end
        x1 = [MatchX(idx,k1) MatchY(idx,k1)];
        x2 = [MatchX(idx,k2) MatchY(idx,k2)];
%         figure;
%         showMatchedFeatures(img{k1},img{k2},x1,x2, 'method', 'montage');
        [~,~,inliners_id] = GetInlinersRANSAC(x1,x2);
%         figure;
%         showMatchedFeatures(img{k1},img{k2},x1,x2, 'method', 'montage');
        Match(setxor(idx,idx(inliners_id)),sum(1:(img_num-1))-sum(1:(img_num-k1))+k2-k1) = 0; % remove outliers
    end
end
save DataRANSAC0.75.mat