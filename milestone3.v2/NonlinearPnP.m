function [C,R] = NonlinearPnP(X,x,K,C0,R0)
% X,x N*3 N*2
q0 = R2q(R0);
params0 = [q0;C0];

opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');
%opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');

fun = @(params)PnPError(X,x,K,params); % user-defined function to compute the vector-valued function
params_optimized = lsqnonlin(fun,params0,[],[],opts);
C = params_optimized(5:7);
R = q2R(params_optimized(1:4));
end

function [error,jacobian] = PnPError(X,x,K,params)
% Measure a PnP error matrix
N = size(X,1); 
P = K*q2R(params(1:4))*[eye(3) -params(5:7)];
error = zeros(2*N,1);

jacobian = zeros(2*N,7);
Rt = q2R(params(1:4));
Ct = params(5:7);
qt = params(1:4);

for i = 1:N
    Xt = X(i,:)';
    error(2*i-1:2*i) = [x(i,1)-(P(1,:)*[Xt;1])/(P(3,:)*[Xt;1]) ; x(i,2)-(P(2,:)*[Xt;1])/(P(3,:)*[Xt;1])];
    
    if nargout > 1
            u = P(1,:)*[Xt;1];
            v = P(2,:)*[Xt;1];
            w = P(3,:)*[Xt;1];
            
            ujacobianC = -[K(1,1)*Rt(1,1)+K(1,3)*Rt(3,1) , K(1,1)*Rt(1,2)+K(1,3)*Rt(3,2) , K(1,1)*Rt(1,3)+K(1,3)*Rt(3,3)];
            vjacobianC = -[K(2,2)*Rt(2,1)+K(2,3)*Rt(3,1) , K(2,2)*Rt(2,2)+K(2,3)*Rt(3,2) , K(2,2)*Rt(2,3)+K(2,3)*Rt(3,3)];
            wjacobianC = - Rt(3,:);
            fjacobianC = [(w*ujacobianC-u*wjacobianC)/(w^2) ; (w*vjacobianC-v*wjacobianC)/(w^2)];
            
            ujacobianR = [K(1,1)*(Xt-Ct)' zeros(1,3) K(1,3)*(Xt-Ct)'];
            vjacobianR = [zeros(1,3) K(2,2)*(Xt-Ct)'  K(2,3)*(Xt-Ct)'];
            wjacobianR = [zeros(1,3) zeros(1,3) (Xt-Ct)'];
            fjacobianR = [(w*ujacobianR-u*wjacobianR)/(w^2) ; (w*vjacobianR-v*wjacobianR)/(w^2)];
            Rjacobianq = [ 0 0 -4*qt(3) -4*qt(4);...
                          -2*qt(4) 2*qt(3) 2*qt(2) -2*qt(1) ;...
                          2*qt(3) 2*qt(4) 2*qt(1) 2*qt(2) ;...
                          2*qt(4) 2*qt(3) 2*qt(2) 2*qt(1) ;...
                           0 -4*qt(2) 0 -4*qt(4) ;...
                          -2*qt(2) -2*qt(1) 2*qt(4) 2*qt(3);...
                          -2*qt(3) 2*qt(4) -2*qt(1) 2*qt(2) ;...
                          2*qt(2) 2*qt(1) 2*qt(4) 2*qt(3) ;...
                           0 -4*qt(2) -4*qt(3) 0];
            fjacobianq = fjacobianR * Rjacobianq;
            
            jacobian(2*i-1:2*i,:)= -[fjacobianq , fjacobianC];
    end
end

end