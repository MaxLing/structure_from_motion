function [im1_pts, im2_pts] = clickCorrespondences(im1, im2) 
% Click on correspondences in two images, and then exit the window
% to return them.
% Inputs:
% im1 - m1 x m1 x 3 image
% im2 - m2 x m2 x 3 image
% Outputs:
% im1_pts - p x 2, where p is the number of clicked points
% im2_pts - p x 2, where p is the number of clicked points

%Ask user to input corresponding points
uiwait(msgbox('Pick corresponding points.'));
[im1_pts, im2_pts] = cpselect(im1, im2,'Wait', true);

%Check if there are equal number of corresponding points
if(size(im1_pts,1)==size(im2_pts,1))
    display('The number of corresponding points are equal');
else
    display('Error!! The number of corresponding points are not equal');
end

end