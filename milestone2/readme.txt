My image, camera intrinsic and correspondence are saved in file /my data.
The provided image, camera intrinsic and correspondence are saved in file /test_images.

Run PlotEpipoles to plot epipolar lines and epipoles.
You need to load camera intrinsic and correspondence mat.file and image first.
Then comment one of two following functions to use another Foundamental matrix.
% F = EstimateFundamentalMatrix(corrs1,corrs2); % for initial estimated F
% [ ~, F, ~, ~, ~ ] = GetTransformation( corrs1, corrs2, K1, K2 ); % for reconstucted F

Run autoCorrespondence to get correspondence before and after RANSAC.
You need to load camera intrinsic mat.file and image first.
You can save corr1,corr2,corr1_new,corr2_new and plot them by runing PlotEpipoles.



Wudao Ling
4/24/2017