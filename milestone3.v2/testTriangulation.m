load('test/matchingRANSAC0.5.mat');
idx = find(Match(:,initialimg1)==1 & Match(:,initialimg2)==1);
x1 = [MatchX(idx,initialimg1) MatchY(idx,initialimg1)];
x2 = [MatchX(idx,initialimg2) MatchY(idx,initialimg2)];

F = EstimateFundamentalMatrix(x1, x2); % Estimate F and E.
E = EssentialMatrixFromFundamentalMatrix(F, K);   
[Cs, Rs] = ExtractCameraPose(E); % Extract C and R.
% Triangulate Points for each C and R.
for i = 1:4
    Xs{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cs{i}, Rs{i}, x1, x2);
end
% Disambiguate the actual pose.
[C, R, X] = DisambiguateCameraPose(Cs, Rs, Xs);
% remove weird points
idx_prime = find(X(:,3) < 0 | abs(X(:,1)) > 100 | abs(X(:,2)) > 100 | abs(X(:,3)) > 100);
X(idx_prime,:) = [];
idx(idx_prime,:) = [];
x1(idx_prime,:) = [];
x2(idx_prime,:) = [];

Cset{1} = zeros(3,1);
Rset{1} = eye(3);
Cset{2} = C;
Rset{2} = R;
PlotCamerasAndPoints( Cset, Rset, X, 1 );

% nonlinear
X = NonlinearTriangulation(K,zeros(3,1),eye(3),C,R,x1,x2,X);
% remove weird points
idx_prime = find(X(:,3) < 0 | abs(X(:,1)) > 100 | abs(X(:,2)) > 100 | abs(X(:,3)) > 100);
X(idx_prime,:) = [];
idx(idx_prime,:) = [];

PlotCamerasAndPoints( Cset, Rset, X, 1 );

idx_frame = [initialimg1,initialimg2]; % index for frame

X3D(idx,:) = X; % idx here is visible 3D point
V(idx, idx_frame) = 1;