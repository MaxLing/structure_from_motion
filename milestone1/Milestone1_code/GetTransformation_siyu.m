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

end
function E = EssentialMatrixFromFundamentalMatrix(F, K1, K2)
E = K2'*F*K1;
[U,D,V] = svd(E);
E = U*diag([1 1 0])*V';
end
function F = EstimateFundamentalMatrix(x1, x2)
A=[x1(:,1).*x2(:,1) x1(:,2).*x2(:,1) x2(:,1) x1(:,1).*x2(:,2) x1(:,2).*x2(:,2) x2(:,2) x1(:,1) x1(:,2), ones(8,1)];

% Get F matrix
[U,D,V] = svd(A);
F=reshape(V(:,9), 3, 3)';
% make rank 2 
[U,D,V] = svd(F);
F=U*diag([D(1,1) D(2,2) 0])*V';
end
function Xset = LinearTriangulation(K1, K2, C1, R1, C2, R2, x1, x2)
P1 = K1*[R1,-R1*C1];
P2 = K2*[R2,-R2*C2];
Xset = zeros(size(x1,1), 3);
for i = 1:size(x1,1)
    skew1 = Vec2Skew([x1(i,:),1]);
    skew2 = Vec2Skew([x2(i,:),1]);
    A = [skew1*P1;skew2*P2];
    [U,D,V] = svd(A);
    X =V(:, end)/V(end,end);
    Xset(i, :) = X(1:3);
end
end
function skew = Vec2Skew(v)
skew = [0, -v(3),v(2);v(3),0,-v(1);-v(2), v(1), 0];
end
function[Cs, Rs] = ExtractCameraPose(E)

[U,D,V] = svd(E);
W = [0,-1,0;
    1,0,0;
    0,0,1];
Rs = cell(4,1);
Cs = cell(4,1);
t1 = U(:,3);
Rs{1} = U*W*V';
Cs{1} = -Rs{1}'*t1;
t2 = -U(:,3);
Rs{2} = U*W*V';
Cs{2} = -Rs{2}'*t2;
t3 = U(:,3);
Rs{3} = U*W'*V';
Cs{3} = -Rs{3}'*t3;
t4 = -U(:,3);
Rs{4} = U*W'*V';
Cs{4} = -Rs{4}'*t4;
for i = 1:4
    if det(Rs{i})<0
        Cs{i} = -Cs{i};
        Rs{i} = -Rs{i};
    end
end

end
function [C, R, X0] = DisambiguateCameraPose(Cs, Rs, Xset)
best_num = 0;
for i = 1:4
    num = 0;
    for j = 1:8
        X = Xset{i}(j,:);
        if Rs{i}(3,:)*(X'-Cs{i})>0
            num = num+1;
        end
        if [0 0 1]*(X')>0
            num = num+1;
        end
    end
    if num > best_num
        C = Cs{i};
        R = Rs{i};
        X0 = Xset{i};
        best_num = num;
    end
end
end