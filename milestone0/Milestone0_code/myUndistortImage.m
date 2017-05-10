function [ I_undist ] = myUndistortImage( I_dist, K, k )
% Given a distorted image, the camera calibration, and distortion
% coefficients, outputs an undistorted version of the image.
%   Inputs:
%   I_dist - (m x n x 3) image matrix
%   K - (3 x 3) intrinsic calibration matrix
%   k - (2 x 1) distortion coefficients
%   Outputs: 
%   I_undist - (m x n x 3) undistorted image matrix
I_dist = im2double(I_dist);
I_undist = zeros(size(I_dist));
[v,u] = find(~isnan(I_undist(:,:,1)));

% u_bar_undistorted = inv(K)* u_undistorted
u_bar_undistorted = K \ [u v ones(length(u),1)]'; % 3*mn
u_bar_distorted = zeros(size(u_bar_undistorted));

% forward distortion
rho2 = u_bar_undistorted(1,:).^2 + u_bar_undistorted(2,:).^2;

u_bar_distorted(1,:) = u_bar_undistorted(1,:).*(1+k(1)*rho2 + k(2)*rho2.^2);
u_bar_distorted(2,:) = u_bar_undistorted(2,:).*(1+k(1)*rho2 + k(2)*rho2.^2);
u_bar_distorted(3,:) = 1;

u_distorted = K*u_bar_distorted; % 3*mn
x_distorted = reshape(u_distorted(1,:),size(I_undist(:,:,1))); % reshape to m*n
y_distorted = reshape(u_distorted(2,:),size(I_undist(:,:,1)));
for i = 1:3
% inverse warping
I_undist(:,:,i) = interp2(I_dist(:,:,i), x_distorted, y_distorted);
end
end