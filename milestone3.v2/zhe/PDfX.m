function PDfX = PDfX(K,R,C,X)
%% Partial differential equation w.r.t World Coordinates in reprojection error problem
X = [X,1];
P = K*R*[eye(3),-C];
x = P*X';%[u;v;w];
PDuX = [K(1,1)*R(1,1)+K(1,3)*R(3,1),K(1,1)*R(1,2)+K(1,3)*R(3,2),K(1,1)*R(1,3)+K(1,3)*R(3,3)];
PDvX = [K(2,2)*R(2,1)+K(2,3)*R(3,2),K(2,2)*R(2,2)+K(2,3)*R(3,2),K(2,2)*R(2,3)+K(2,3)*R(3,3)];
PDwX = [R(3,1),R(3,2),R(3,3)];
PDfX = [(x(3)*PDuX-x(1)*PDwX)/x(3)^2;
        (x(3)*PDvX-x(2)*PDwX)/x(3)^2];
end