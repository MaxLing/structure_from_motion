function [X] = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)

%   Inputs:
%   x1: (8x2) points in Image 1 N*2
%   x2: (8x2) corresponding points in Image 2
%   K1: (3x2) intrinsic calibration for Camera 1
%   K2: (3x2) intrinsic calibration for Camera 2
%   C1, R1 first camera pose
%   C2, R2 second camera pose
%   Output: N*3 matrix representing 3D triangulated point

N = size(x1,1); % num of correspondences
X = zeros(N,3);

P1 = K*R1*[eye(3), -C1];
P2 = K*R2*[eye(3), -C2];

for i = 1:N
    X1=[x1(i,:) ,1]; 
    X2=[x2(i,:) ,1];  
    S1 =[0 -X1(3) X1(2); X1(3) 0 -X1(1); -X1(2) X1(1) 0];
    S2 =[0 -X2(3) X2(2); X2(3) 0 -X2(1); -X2(2) X2(1) 0];
    F=[S1*P1;S2*P2];
    [U,D,V] = svd(F);
    X_prime = (V(:,end)/V(end,end))';
    X(i,:) = X_prime(1:3);
end

end