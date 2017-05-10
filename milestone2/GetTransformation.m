function [ E, F, C, R, X0 ] = GetTransformation( x1, x2, K1, K2 )
%GETTRANSFORMATION Computes the transformation between the two cameraman.
%   Inputs:
%   x1: (8x2) points in Image 1
%   x2: (8x2) corresponding points in Image 2
%   K1: (3x2) intrinsic calibration for Camera 1
%   K2: (3x2) intrinsic calibration for Camera 2
%   Outputs:
%   E: (3x3) essential matrix
%   F: (3x3) fundamental matrix
%   C: (3x1) transformation from Camera 1 to Camera 2
%   R: (3x3) rotation from Camera 1 to Camera 2
%   X0: (nx3) the correspondence points in 3D

% Estimate F and E.
F = EstimateFundamentalMatrix(x1, x2);
E = EssentialMatrixFromFundamentalMatrix(F, K1, K2);   

% Extract C and R.
[Cs, Rs] = ExtractCameraPose(E);

%Triangulation
fprintf('Triangulating points between first two images...');
fprintf('\n');

% Triangulate Points for each C and R.
for i = 1:4
    Xset{i} = LinearTriangulation(K1, K2, zeros(3,1), eye(3), Cs{i}, Rs{i}, x1, x2);
end

% Disambiguate the actual pose.
[C, R, X0] = DisambiguateCameraPose(Cs, Rs, Xset);

% reconstruct F
t = -R*C;
F = (K2'^-1)*[0, -t(3),t(2);t(3),0,-t(1);-t(2), t(1), 0]*R*K1^-1;

end

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

function [E] = EssentialMatrixFromFundamentalMatrix(F,K1,K2)

% Input: F K1 K2
% Output: E

E = K2'*F*K1; % 

% 2nd cleanup
[U,D,V]=svd(E);
D = [1 0 0;0 1 0;0 0 0];
E = U*D*V';

end

function [Cs, Rs] = ExtractCameraPose(E)
% Input: E
% Output: 4 camera pose configurations

Rs = cell(4,1);
Cs = cell(4,1);

[U,D,V]=svd(E);

% 2 set of t
t1 = U(:,3);
t2 = -U(:,3);

% 2 set of W
W1 = [0,-1,0;1 0 0;0 0 1];
W2 = W1';

Rs{1} = U*W1*V';
Cs{1} = -Rs{1}'*t1;
Rs{2} = U*W1*V';
Cs{2} = -Rs{2}'*t2;
Rs{3} = U*W2*V';
Cs{3} = -Rs{3}'*t1;
Rs{4} = U*W2*V';
Cs{4} = -Rs{4}'*t2;
for i = 1:4
    if det(Rs{i})<0
        Rs{i} = -Rs{i};
        Cs{i} = -Cs{i};
    end
end


end

function [X] = LinearTriangulation(K1, K2, C1, R1, C2, R2, x1, x2)

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

P1 = K1*R1*[eye(3), -C1];
P2 = K2*R2*[eye(3), -C2];

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

function [C, R, X0] = DisambiguateCameraPose(Cs, Rs, Xset)
% input: 
% Cs,Rs: four camera configurations
% Xset: four sets of triangulation points
% output:
% C,R: the correct camera pose
% X0: the correct triangulation points
max_cheirality = 0;
best_id = NaN;

for i=1:4
    C = Cs{i};
    r3 = Rs{i}(3,:);
    Xs = Xset{i};
    cheirality =0;
    for j = 1:8
        if r3*(Xs(j,:)'-C) > 0 % triangulation should be in front of the camera 2
            cheirality=cheirality+1;
        end
        if [0 0 1]*(Xs(j,:)') > 0 % triangulation should be in front of the camera 1
            cheirality=cheirality+1;
        end
    end
    if cheirality > max_cheirality
        max_cheirality = cheirality;
        best_id = i;
    end
end
C = Cs{best_id};
R = Rs{best_id};
X0 = Xset{best_id};

end

