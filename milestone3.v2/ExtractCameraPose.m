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