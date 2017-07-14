function [F] = EstimateFundamentalMatrix(x1,x2)

% Input: x1 x2 8*2 matrix that repersent correspondences each row
% Output: F
N = size(x1,1);

A = [x1(:,1).*x2(:,1) x1(:,2).*x2(:,1) x2(:,1) x1(:,1).*x2(:,2) x1(:,2).*x2(:,2) x2(:,2) x1(:,1) x1(:,2) ones(N,1)];
[~,~,V]=svd(A);
f = V(:,end);
F = [f(1:3)';f(4:6)';f(7:9)'];

% 1st cleanup rank = 2
[u,d,v] = svd(F);
d(3,3)=0;
F = u*d*v';
end