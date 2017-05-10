function [lines] = myEpipolarLine(F,points)

% input F 3*3 foundamental matrix map from camera B to camera A
% input points n*2 points in camera B
% output lines n*3 epipolar lines in camera A corresponding to points in B

n = size(points,1); 
points = [points,ones(n,1)]'; % 3*n
lines = zeros(n,3);

for i = 1: n
    lines(i,:) = (F*points(:,i))';
end

end