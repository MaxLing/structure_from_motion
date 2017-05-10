function [c1,c2,best_inliners_id] = GetInlinersRANSAC(corr1,corr2)

% corr1,corr2: n corner points before/after RANSAC
% num: num of points needed for random sampling
% rate: success rate of RANSAC, which decide the iteration
% dist: threshold dist between inliners and fitting line
% ratio: threshold of ratio, inliners to all points
num = 8;
rate = 0.99;
dist = 0.75; %0.5 0.8
ratio = 0.4;

m = size(corr1,1); % num of all points
iter = round(log(1-rate)/log(1-ratio^num)); % too small here???
max_inliners_num = 0; 

for i = 1:iter
    % randomly select num points for sampling
    random = randperm(m,num); 
    x1 = corr1(random,:);
    x2 = corr2(random,:);
    F = EstimateFundamentalMatrix(x1,x2);
    
    % threshold the inliers
    
    distance = zeros(m,1);
    for j = 1:m
%          distance(j) = abs([corr2(j,:),1]*F*[corr1(j,:),1]');
        distance(j) = abs([corr2(j,:),1]*F*[corr1(j,:),1]')/norm(F(1:2,:)*[corr1(j,:),1]');
        % dist = abs(V^T*F*U)/sqrt((F1*U)^2+(F2*U)^2)
    end
    inliners_id = find(distance<=dist);
    inliners_num = length(inliners_id);
    
    % Update inliers and the model if a better model is found     
    if inliners_num > max_inliners_num
         max_inliners_num = inliners_num;
         best_inliners_id = inliners_id;
    end
end

c1 = corr1(best_inliners_id,:);
c2 = corr2(best_inliners_id,:);
   
end