num = 6;
rate = 0.99;
dist = 10;
ratio = 0.5;

m = size(X,1); % num of all points
iter = round(log(1-rate)/log(1-ratio^num)); % too small here???
max_inliners_num = 0; 
best_dist = [];

for i = 1:iter
    % randomly select num points for sampling
    random = randperm(m,num); 
    xs = x(random,:);
    Xs = X(random,:);
    [Cs,Rs]=LinearPnP(Xs,xs,K);
    Ps = K*Rs*[eye(3) -Cs];
    
    % threshold the inliers
    distance = zeros(m,1);
    for j = 1:m
        distance(j) = (x(j,1)-(Ps(1,:)*[X(j,:)';1])/(Ps(3,:)*[X(j,:)';1]))^2 + (x(j,2)-(Ps(2,:)*[X(j,:)';1])/(Ps(3,:)*[X(j,:)';1]))^2;
    end
    inliners_id = find(distance<=dist);
    inliners_num = length(inliners_id);
    
    % Update inliers and the model if a better model is found     
    if inliners_num > max_inliners_num
         max_inliners_num = inliners_num;
         best_inliners_id = inliners_id;
         best_dist = distance;
    end
end

Xbest = X(best_inliners_id,:);
xbest = x(best_inliners_id,:);
[C,R]=LinearPnP(Xbest,xbest,K);