function PDfC = PDfC(K,R,C,X)
%% Partial differential equation w.r.t Camera Center in reprojection error problem
X = [X,1];
P = K*R*[eye(3),-C];
x = P*X';%[u;v;w];
PDuC = -[K(1,1)*R(1,1)+K(1,3)*R(3,1),K(1,1)*R(1,2)+K(1,3)*R(3,2),K(1,1)*R(1,3)+K(1,3)*R(3,3)];
PDvC = -[K(2,2)*R(2,1)+K(2,3)*R(3,2),K(2,2)*R(2,2)+K(2,3)*R(3,2),K(2,2)*R(2,3)+K(2,3)*R(3,3)];
PDwC = -[R(3,1),R(3,2),R(3,3)];
PDfC = [(x(3)*PDuC-x(1)*PDwC)/x(3)^2;
        (x(3)*PDvC-x(2)*PDwC)/x(3)^2];
end