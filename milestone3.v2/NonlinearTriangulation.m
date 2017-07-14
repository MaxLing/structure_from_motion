function X = NonlinearTriangulation(K,C1,R1,C2,R2,x1,x2,X0)

% C1 R1 first camera pose
% C2 R2 second camera pose
% x1 x2 N*2 correspondence
% X X0 N*3 3d point

P1 = K*R1*[eye(3) -C1];
P2 = K*R2*[eye(3) -C2];
X = zeros(size(X0,1),3);
opts = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');
%opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true, 'Algorithm', 'levenberg-marquardt', 'TolX', 1e-64, 'TolFun', 1e-64, 'MaxFunEvals', 1e64, 'MaxIter', 1e64, 'Display', 'iter');

for i =1:size(X0,1)
params0 = X0(i,:)'; % 3*1
fun = @(params) TriangulationError(P1,P2,R1,R2,K,x1(i,:),x2(i,:),params); % user-defined function to compute the vector-valued function
X(i,:) = (lsqnonlin(fun,params0,[],[],opts))';
end

end

function [error,jacobian] = TriangulationError(P1,P2,R1,R2,K,x1,x2,X) 
% Measure a triangulation error matrix
    error = [x1(1)-(P1(1,:)*[X;1])/(P1(3,:)*[X;1]);...
             x1(2)-(P1(2,:)*[X;1])/(P1(3,:)*[X;1]);... % camera 1 
             x2(1)-(P2(1,:)*[X;1])/(P2(3,:)*[X;1]);...
             x2(2)-(P2(2,:)*[X;1])/(P2(3,:)*[X;1])];   % camera 2
    
    if nargout > 1
        u1 = P1(1,:)*[X;1];
        v1 = P1(2,:)*[X;1];
        w1 = P1(3,:)*[X;1];
        ujacobianX1 = [K(1,1)*R1(1,1)+K(1,3)*R1(3,1) , K(1,1)*R1(1,2)+K(1,3)*R1(3,2) , K(1,1)*R1(1,3)+K(1,3)*R1(3,3)];
        vjacobianX1 = [K(2,2)*R1(2,1)+K(2,3)*R1(3,1) , K(2,2)*R1(2,2)+K(2,3)*R1(3,2) , K(2,2)*R1(2,3)+K(2,3)*R1(3,3)];
        wjacobianX1 =  R1(3,:);
        fjacobianX1 = [(w1*ujacobianX1-u1*wjacobianX1)/(w1^2) ; (w1*vjacobianX1-v1*wjacobianX1)/(w1^2)];
        
        u2 = P2(1,:)*[X;1];
        v2 = P2(2,:)*[X;1];
        w2 = P2(3,:)*[X;1];
        ujacobianX2 = [K(1,1)*R2(1,1)+K(1,3)*R2(3,1) , K(1,1)*R2(1,2)+K(1,3)*R2(3,2) , K(1,1)*R2(1,3)+K(1,3)*R2(3,3)];
        vjacobianX2 = [K(2,2)*R2(2,1)+K(2,3)*R2(3,1) , K(2,2)*R2(2,2)+K(2,3)*R2(3,2) , K(2,2)*R2(2,3)+K(2,3)*R2(3,3)];
        wjacobianX2 =  R2(3,:);
        fjacobianX2 = [(w2*ujacobianX2-u2*wjacobianX2)/(w2^2) ; (w2*vjacobianX2-v2*wjacobianX2)/(w2^2)];
        
        jacobian = -[fjacobianX1;fjacobianX2];
    end
end