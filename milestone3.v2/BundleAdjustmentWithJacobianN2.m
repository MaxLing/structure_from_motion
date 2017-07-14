 function [Cset, Rset, Xset] = BundleAdjustmentWithJacobianN2(Cset,Rset,X,K,x,y,V)
% BundleAdjustmentWithJacobian
% without MATLAB built-in nonlinear optimization
% By Wudao Ling 6/13/2017

nIters = 5;
lambda = 10000000;
I = length(Cset); % num of cameras
J = size(X,1); % num of points

% initialization
U = sum(reshape(V,[],1)); % num of visible points
params = zeros(7*(I-1)+3*J,1);
params(7*(I-1)+1:7*(I-1)+3*J) = reshape(X',3*J,1);
for i = 2:I % we don't tune camera 1
    q = R2q(Rset{i});
    params((i-2)*7+1:(i-2)*7+4) = q;
    params((i-2)*7+5:(i-2)*7+7) = Cset{i};
end

% for i = 1:I
%         if i==1
%             Rtest = [1 0 0; 0 1 0; 0 0 1];
%             Ctest = [0;0;0];
%         else
%             Rtest = q2R(params((i-2)*7+1:(i-2)*7+4));
%             Ctest = params((i-2)*7+5:(i-2)*7+7);
%         end
%     Xtest = (reshape(params(7*(I-1)+1:7*(I-1)+3*J),3,[]))' ;
%     gerror = TestReprojection(Ctest,Rtest,Xtest,K,x,y,V,i);
%     GerrorProcess(i,1)=gerror;
% end

% update
for k = 1:nIters
    [error,Jacobian] = BundleError(K,x,y,V,I,J,U,params);
    delta = (Jacobian'*Jacobian+lambda*eye(size(Jacobian,2)))\(Jacobian')*error;
    params = params + delta;
    for i = 2:I
        params((i-2)*7+1:(i-2)*7+4) = params((i-2)*7+1:(i-2)*7+4)/norm(params((i-2)*7+1:(i-2)*7+4)); % normalize quaternion
    end
    % test cams
%     for i = 1:I
%         if i==1
%             Rtest = [1 0 0; 0 1 0; 0 0 1];
%             Ctest = [0;0;0];
%         else
%             Rtest = q2R(params((i-2)*7+1:(i-2)*7+4));
%             Ctest = params((i-2)*7+5:(i-2)*7+7);
%         end
%     Xtest = (reshape(params(7*(I-1)+1:7*(I-1)+3*J),3,[]))' ;
%     gerror = TestReprojection(Ctest,Rtest,Xtest,K,x,y,V,i);
%     GerrorProcess(i,k+1)=gerror;
%     end
end

% output
Cset = cell([1 I]);
Rset = cell([1 I]);
Rset{1} = [1 0 0; 0 1 0; 0 0 1];
Cset{1} = [0;0;0];
for i=2:I
    Rset{i} = q2R(params((i-2)*7+1:(i-2)*7+4));
    Cset{i} = params((i-2)*7+5:(i-2)*7+7);
end
Xset = (reshape(params(7*(I-1)+1:7*(I-1)+3*J),3,[]))';
end

function [error,jacobian] = BundleError(K,x,y,V,I,J,U,params)
% build a error and jacobian matrix for bundle adjustment
error = zeros(U*2,1);
jacobian = zeros(U*2,7*(I-1)+3*J);
k = 0;

for j = 1:J
    Xt = params(7*(I-1)+(j-1)*3+1:7*(I-1)+(j-1)*3+3);
    for i = 1:I
        if (i==1)
            qt = [1;0;0;0];
            Ct = [0;0;0];
        else
            qt = params((i-2)*7+1:(i-2)*7+4);
            Ct = params((i-2)*7+5:(i-2)*7+7);
        end
        Rt = q2R(qt);
        p = K*Rt*[eye(3),-Ct];
        if (V(j,i)==1)
            k = k+1;
            error((k-1)*2+1:(k-1)*2+2) = [x(j,i) - (p(1,:)*[Xt;1])/(p(3,:)*[Xt;1]);...
                                          y(j,i) - (p(2,:)*[Xt;1])/(p(3,:)*[Xt;1])];                                      
            u = p(1,:)*[Xt;1];
            v = p(2,:)*[Xt;1];
            w = p(3,:)*[Xt;1];
            
            ujacobianX = [K(1,1)*Rt(1,1)+K(1,3)*Rt(3,1) , K(1,1)*Rt(1,2)+K(1,3)*Rt(3,2) , K(1,1)*Rt(1,3)+K(1,3)*Rt(3,3)];
            vjacobianX = [K(2,2)*Rt(2,1)+K(2,3)*Rt(3,1) , K(2,2)*Rt(2,2)+K(2,3)*Rt(3,2) , K(2,2)*Rt(2,3)+K(2,3)*Rt(3,3)];
            wjacobianX = Rt(3,:);
            fjacobianX = [(w*ujacobianX-u*wjacobianX)/(w^2) ; (w*vjacobianX-v*wjacobianX)/(w^2)];
            
            jacobian((k-1)*2+1:(k-1)*2+2,7*(I-1)+(j-1)*3+1:7*(I-1)+(j-1)*3+3)= fjacobianX;
            
            if (i~=1)
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
            
            jacobian((k-1)*2+1:(k-1)*2+2,(i-2)*7+1:(i-2)*7+7)= [fjacobianq , fjacobianC];
            end
        end
    end
end
jacobian = sparse(jacobian);
end