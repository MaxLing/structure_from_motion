function [Cset, Rset, Xset] = BundleAdjustmentWithJacobian(Cset,Rset,X,K,x,y,V)

I = length(Cset);
J = size(X,1);

params0 = zeros(7*I+3*J,1);
for i = 1:I
    q = R2q(Rset{i});
    params0((i-1)*7+1) = q(1);
    params0((i-1)*7+2) = q(2);
    params0((i-1)*7+3) = q(3);
    params0((i-1)*7+4) = q(4);
    params0((i-1)*7+5) = Cset{i}(1);
    params0((i-1)*7+6) = Cset{i}(2);
    params0((i-1)*7+7) = Cset{i}(3);
end
for j=1:J
    params0(I*7+(j-1)*3+1) = X(j,1);
    params0(I*7+(j-1)*3+2) = X(j,2);
    params0(I*7+(j-1)*3+3) = X(j,3);
end


%opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt','TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'Display', 'iter');
opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true,'Algorithm', 'levenberg-marquardt','TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'Display', 'iter');

fun = @(params)BundleError(K,x,y,V,I,J,params); 
params_optimized = lsqnonlin(fun,params0,[],[],opts);

[Cset, Rset, Xset] = decodeParams(params_optimized,I,J);
end

function [error,jacobian] = BundleError(K,x,y,V,I,J,params)
% build a error and jacobian matrix for bundle
[Cset, Rset, X] = decodeParams(params, I, J);
error = zeros(I*J*2,1);
jacobian = zeros(I*J*2,7*I+3*J);
for i = 1:I
    P = K*Rset{i}*[eye(3),-Cset{i}];
    for j = 1:J
        Xt = X(j,:)';
        error(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2) = [V(j,i)*(x(j,i) - (P(1,:)*[Xt;1])/(P(3,:)*[Xt;1]));...
                                                      V(j,i)*(y(j,i) - (P(2,:)*[Xt;1])/(P(3,:)*[Xt;1]))];
        if nargout > 1
            Rt = Rset{i};
            Ct = Cset{i};
            qt = R2q(Rt);
            u = P(1,:)*[Xt;1];
            v = P(2,:)*[Xt;1];
            w = P(3,:)*[Xt;1];
            
            ujacobianC = -[K(1,1)*Rt(1,1)+K(1,3)*Rt(3,1) , K(1,1)*Rt(1,2)+K(1,3)*Rt(3,2) , K(1,1)*Rt(1,3)+K(1,3)*Rt(3,3)];
            vjacobianC = -[K(2,2)*Rt(2,1)+K(2,3)*Rt(3,1) , K(2,2)*Rt(2,2)+K(2,3)*Rt(3,2) , K(2,2)*Rt(2,3)+K(2,3)*Rt(3,3)];
            wjacobianC = -Rt(3,:);
            fjacobianC = [V(j,i)*(w*ujacobianC-u*wjacobianC)/(w^2) ; V(j,i)*(w*vjacobianC-v*wjacobianC)/(w^2)];

            ujacobianX = -ujacobianC;
            vjacobianX = -vjacobianC;
            wjacobianX = -wjacobianC;
            fjacobianX = [V(j,i)*(w*ujacobianX-u*wjacobianX)/(w^2) ; V(j,i)*(w*vjacobianX-v*wjacobianX)/(w^2)];
            
            ujacobianR = [K(1,1)*(Xt-Ct)' zeros(1,3) K(1,3)*(Xt-Ct)'];
            vjacobianR = [zeros(1,3) K(2,2)*(Xt-Ct)'  K(2,3)*(Xt-Ct)'];
            wjacobianR = [zeros(1,3) zeros(1,3) (Xt-Ct)'];
            fjacobianR = [V(j,i)*(w*ujacobianR-u*wjacobianR)/(w^2) ; V(j,i)*(w*vjacobianR-v*wjacobianR)/(w^2)];
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
            
            jacobian(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2,(i-1)*7+1:(i-1)*7+7)= -[fjacobianq , fjacobianC];
            jacobian(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2,7*I+(j-1)*3+1:7*I+(j-1)*3+3)= -fjacobianX;
            
        end
    end
end
jacobian = sparse(jacobian);
end


function [Cset, Rset, X] = decodeParams(p, I, J)

Cset = cell([1 I]);
Rset = cell([1 I]);
X = zeros(J,3);

for i=1:I
    Rset{i} = q2R(p((i-1)*7+1:(i-1)*7+4));
    Cset{i} = p((i-1)*7+5:(i-1)*7+7);
end

for j=1:J
    X(j,:) = (p(I*7+(j-1)*3+1:I*7+(j-1)*3+3))';
end
end