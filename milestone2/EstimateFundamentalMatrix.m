function [F] = EstimateFundamentalMatrix(x1,x2,K1,K2)

% Input: x1 x2 n*2 matrix that repersent correspondences each row
% Output: F
N = size(x1,1);

A = [x1(:,1).*x2(:,1) x1(:,2).*x2(:,1) x2(:,1) x1(:,1).*x2(:,2) x1(:,2).*x2(:,2) x2(:,2) x1(:,1) x1(:,2) ones(N,1)];
[~,~,V]=svd(A);
F = reshape(V(:,end),3,3)';

% 1st cleanup rank = 2
[u,d,v] = svd(F);
d(3,3)=0;
F = u*d*v';

% E = K2'*F*K1; 
% % 2nd cleanup rank = 2
% [U,D,V]=svd(E);
% D = [1 0 0;0 1 0;0 0 0]; % problem here!!!!!
% % D(3,3)=0;
% E = U*D*V';
% 
% F = inv(K2')*E*inv(K1);

end