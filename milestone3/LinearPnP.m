function [C,R]=LinearPnP(X,x,K)

% X,x N*3 N*2
N = size(X,1);
x = (K \ [x, ones(N,1)]')'; % N*3 normalized x
X = [X, ones(N,1)]; % N*4

% solve PnP
A = [X, zeros(N,4), -(x(:,1)*[1 1 1 1]).*X;...
     zeros(N,4), X, -(x(:,2)*[1 1 1 1]).*X];

 [U,D,V]=svd(A);
 P = V(:,end);
 P = [P(1:4)';P(5:8)';P(9:12)'];
%  P = K \ P; 
 
 % decomposition
 [u,d,v] = svd(P(:,1:3));
 R = u*v';
 t = P(:,4)/d(1,1);
 C = -R'*t;
 
if det(R)<0
        R = -R;
        C = -C;
end

end