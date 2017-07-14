function [C, R, X0] = DisambiguateCameraPose(Cs, Rs, Xset)
% input: 
% Cs,Rs: four camera configurations
% Xset: four sets of triangulation points
% output:
% C,R: the correct camera pose
% X0: the correct triangulation points
max_cheirality = 0;
best_id = NaN;

for i=1:4
    C = Cs{i};
    r3 = Rs{i}(3,:);
    Xs = Xset{i};
    cheirality =0;
    for j = 1:8
        if r3*(Xs(j,:)'-C) > 0 % triangulation should be in front of the camera 2
            cheirality=cheirality+1;
        end
        if [0 0 1]*(Xs(j,:)') > 0 % triangulation should be in front of the camera 1
            cheirality=cheirality+1;
        end
    end
    if cheirality > max_cheirality
        max_cheirality = cheirality;
        best_id = i;
    end
end
C = Cs{best_id};
R = Rs{best_id};
X0 = Xset{best_id};

end