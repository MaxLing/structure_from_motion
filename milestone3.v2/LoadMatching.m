function [Mx,My,M] = LoadMatching(filename, image_idx, img_num)

fid = fopen(filename);
fscanf(fid, '%s', 1);
n = fscanf(fid, '%d', 1); % nFeatures
Mx = zeros(n, img_num);
My = zeros(n, img_num);
M = zeros(n, sum(1:img_num-1));
for i = 1 : n
    m = fscanf(fid, '%d', 1); % all matches for ith feature
    fscanf(fid, '%d', 3); % drop RGB
    Mx(i, image_idx) = fscanf(fid, '%f', 1); % u
    My(i, image_idx) = fscanf(fid, '%f', 1); % v
%     M(i, image_idx) = 1; % current image
    for k = 1 : m-1
        j = fscanf(fid, '%d', 1); % image id
        Mx(i, j) = fscanf(fid, '%f', 1); % u
        My(i, j) = fscanf(fid, '%f', 1); % v
%     M(i, j) = 1; % following image
        M(i, sum(1:(img_num-1))-sum(1:(img_num-image_idx))+j-image_idx) = 1; % current image + following images pairs
    end
end
fclose(fid);

end