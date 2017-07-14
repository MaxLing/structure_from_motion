function [Cset, Rset, Xset] = BundleAdjustmentWithJacobianN(Cset,Rset,X,K,x,y,V)
% BundleAdjustmentWithJacobian
% without MATLAB built-in nonlinear optimization
% By Wudao Ling 6/12/2017

nIters = 10;
I = length(Cset); % num of cameras
J = size(X,1); % num of points

% initialization
P = zeros(7*(I-1),1);
X = reshape(X',3*J,1);
for i = 2:I % we don't tune camera 1
    q = R2q(Rset{i});
    P((i-2)*7+1:(i-2)*7+4) = q;
    P((i-2)*7+5:(i-2)*7+7) = Cset{i};
end

% update
for k = 1:nIters
    [error,JacobianP,JacobianX] = BundleError(K,x,y,V,I,J,P,X);
    A = JacobianP'*JacobianP;
    B = JacobianP'*JacobianX;
    D = JacobianX'*JacobianX;
    eP = JacobianP'*error;
    eX = JacobianX'*error;
    delta_P = (A-B/D*B')\(eP-B/D*eX);
    delta_X = D\(eX-B'*delta_P); 
%     delta_p = inv(A-B*inv(D)*B')*(eP-B*inv(D)*eX);
%     delta_X = inv(D)*(eX-B'*delta_P); 
    X = X + delta_X;
    P = P + delta_P;
    for i = 2:I
        P((i-2)*7+1:(i-2)*7+4) = P((i-2)*7+1:(i-2)*7+4)/norm(P((i-2)*7+1:(i-2)*7+4)); % normalize quaternion
    end
end

% output
Cset = cell([1 I]);
Rset = cell([1 I]);
Xset = zeros(J,3);
Rset{1} = [1 0 0; 0 1 0; 0 0 1];
Cset{1} = [0;0;0];
for i=2:I
    Rset{i} = q2R(P((i-2)*7+1:(i-2)*7+4));
    Cset{i} = P((i-2)*7+5:(i-2)*7+7);
end
for j=1:J
    Xset(j,:) = (X((j-1)*3+1:(j-1)*3+3))';
end
end

function [error,jacobianP,jacobianX] = BundleError(K,x,y,V,I,J,P,X)
% build a error and 2 jacobian matrix for bundle adjustment
error = zeros(I*J*2,1);
jacobianP = zeros(I*J*2,7*(I-1));
jacobianX = zeros(I*J*2,3*J);
for i = 1:I
    if (i==1)
        qt = [1;0;0;0];
        Ct = [0;0;0];
    else
        qt = P((i-2)*7+1:(i-2)*7+4);
        Ct = P((i-2)*7+5:(i-2)*7+7);
    end
    Rt = q2R(qt);
    p = K*Rt*[eye(3),-Ct];
    
    for j = 1:J
        Xt = X((j-1)*3+1:(j-1)*3+3);
        if (V(j,i)==1)
            error(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2) = [x(j,i) - (p(1,:)*[Xt;1])/(p(3,:)*[Xt;1]);...
                                                          y(j,i) - (p(2,:)*[Xt;1])/(p(3,:)*[Xt;1])];                                      
            u = p(1,:)*[Xt;1];
            v = p(2,:)*[Xt;1];
            w = p(3,:)*[Xt;1];
            
            ujacobianX = [K(1,1)*Rt(1,1)+K(1,3)*Rt(3,1) , K(1,1)*Rt(1,2)+K(1,3)*Rt(3,2) , K(1,1)*Rt(1,3)+K(1,3)*Rt(3,3)];
            vjacobianX = [K(2,2)*Rt(2,1)+K(2,3)*Rt(3,1) , K(2,2)*Rt(2,2)+K(2,3)*Rt(3,2) , K(2,2)*Rt(2,3)+K(2,3)*Rt(3,3)];
            wjacobianX = Rt(3,:);
            fjacobianX = [(w*ujacobianX-u*wjacobianX)/(w^2) ; (w*vjacobianX-v*wjacobianX)/(w^2)];
            
%             ujacobianC = -ujacobianX;
%             vjacobianC = -vjacobianX;
%             wjacobianC = -wjacobianX;
%             fjacobianC = [(w*ujacobianC-u*wjacobianC)/(w^2) ;(w*vjacobianC-v*wjacobianC)/(w^2)];
            fjacobianC = -fjacobianX;
            
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
            
            jacobianX(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2,(j-1)*3+1:(j-1)*3+3)= fjacobianX;
            if (i~=1)
                jacobianP(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2,(i-2)*7+1:(i-2)*7+7)= [fjacobianq , fjacobianC];
            end
        end
    end
end
jacobianP = sparse(jacobianP);
jacobianX = sparse(jacobianX);
end