function PDfq = PDfq(K,R,C,X)
%% Partial differential equation w.r.t quaternion in reprojection error problem
X = [X,1];
P = K*R*[eye(3),-C];
x = P*X';%[u;v;w];
X(:,4) = [];
PDuR = [K(1,1)*(X-C'),zeros(1,3),K(1,3)*(X-C')];
PDvR = [zeros(1,3),K(2,2)*(X-C'),K(2,3)*(X-C')];
PDwR = [zeros(1,6),(X-C')];
PDfR = [(x(3)*PDuR-x(1)*PDwR)/x(3)^2;
        (x(3)*PDvR-x(2)*PDwR)/x(3)^2];
q = R2q(R);
PDRq = [0,0,-4*q(3),-4*q(4);
        2*q(4),2*q(3),2*q(2),-2*q(1);
        2*q(3),2*q(4),2*q(1),2*q(2);
        2*q(4),2*q(3),2*q(2),2*q(1);
        0,-4*q(2),0,-4*q(4);
        -2*q(2),-2*q(1),2*q(4),2*q(3);
        -2*q(3),2*q(4),-2*q(1),2*q(2);
        2*q(2),2*q(1),2*q(4),2*q(3);
        0,-4*q(2),-4*q(3),0];
PDfq = PDfR*PDRq;
end