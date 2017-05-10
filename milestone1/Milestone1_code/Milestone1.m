%% TODO: You should fill in these values appropriately.
% You should run this script for each pair of images you take.

% It might be easiest if you save these in a .mat during Milestone 0,
% then load them here.
load('correspondence_ransac0.005.mat') % change here for different pairs
% corresondences1.mat for case1,corresondences2.mat for case2,corresondences3.mat for case3
correspondence_points_1 = corrs1; % Should be p x 2, where p is number of points
correspondence_points_2 = corrs2; % Should be p x 2, where p is number of points

K1 = [2.768201922155784e+03,-3.283946709730422,1.638682157182430e+03;0,2.767657671862577e+03,1.244945553677645e+03;0,0,1]; 
% Should be 3 x 3, calibration matrix for first camera
K2 = [3.456597718123323e+03, 4.948208007856249,2.009097474856611e+03;0,3.465840808785538e+03,1.535384201592209e+03;0,0,1]; 
% Should be 3 x 3, calibration matrix for second camera

k1 = [0.083541099963479;-0.453999314312728]; % Should be 2 x 1, radial distortion for first camera
k2 = [0.036361182076967; 0.114666012903610]; % Should be 2 x 1, radial distortion for second camera

% from test images
% K = [568.996140852000,0,643.210559410000;0,568.988362396000,477.982801038000;0,0,1];
% K1=K;
% K2=K;
% k1=[0;0];
% k2=[0;0];

%% Do the computation
[ E, F, C, R, X ] = GetTransformation( ...
    correspondence_points_1, correspondence_points_2, K1, K2 );

% You'll want to modify this to save a separate .mat for each pair of images.
save transform.mat E F C R X

%% Visualization

% Assume the first cameraman is at the origin
Cset{1} = [0;0;0];
Rset{1} = eye(3);

Cset{2} = C;
Rset{2} = R;

% You should include the output figure in your report for each pair of
% images. You may need to adjust the scale depending on how are are your
% cameramen are from each other.
PlotCamerasAndPoints(Cset, Rset, X, 25);