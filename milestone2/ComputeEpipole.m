function [e1, e2] = ComputeEpipole(F)

e1 = null(F);
e2 = null(F');
e1 = e1/e1(3);
e2 = e2/e2(3);

% [~,~,V1] = svd(F);
% e1 = V1(:,end); % 3*1
% 
% [~,~,V2] = svd(F');
% e2 = V2(:,end);
end