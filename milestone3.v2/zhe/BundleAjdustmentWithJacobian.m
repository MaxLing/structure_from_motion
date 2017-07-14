function [Cset,Rset,Xset] = BundleAjdustmentWithJacobian(Cset,Rset,X,K,traj,V)
%% Bundle Adjustment

n = size(Cset,2);
flag = [];
X_init = [];
for i = 1:n
    if ~isempty(Rset{i})
        q = R2q(Rset{i});
        X_init = [X_init;q;Cset{i}];
        flag = [flag;i];
    end
end
X0 = X';
X_init = [X_init;X0(:)];

Cam = size(flag,1);

opts = optimoptions(@lsqnonlin,'Jacobian','on','display','iter');
fun = @(params)reproError(K,traj,V,params,Cam);
X_opt = lsqnonlin(fun,X_init,[],[],opts);

m = size(flag,1);
for i = 1:m
    R_prime = q2R(X_opt(7*i-6:7*i-3)/norm(X_opt(7*i-6:7*i-3)));
    C_prime = X_opt(7*i-2:7*i);
    ind = flag(i,1);
    Rset{ind} = R_prime;
    Cset{ind} = C_prime;
end
Xset = X_opt(m*7+1:end);
n = size(Xset,1);
Xset = reshape(Xset,[3,n/3]);
Xset = Xset';
end

function [e,J] = reproError(K,traj,V,params,Cam)
e = [];
X0 = params(7*Cam+1:end);
n = size(X0,1);
X0 = reshape(X0,[3,n/3]);
X0 = X0';
n = size(X0,1);
X0 = [X0,ones(n,1)];
J = sparse(2*Cam*n,7*Cam+3*n);
for i = 1:Cam
    q = params(7*i-6:7*i-3);
    R = q2R(q/norm(q));
    C = params(7*i-2:7*i);
    P = K*R*[eye(3),-C];
    x_repro = P*X0';%[u;v;w];
    x = [];
    for j = 1:n
        Jc = [PDfq(K,R,C,X0(j,1:3)),PDfC(K,R,C,X0(j,1:3))];
        Jc = Jc';
        Jp = PDfX(K,R,C,X0(j,1:3));
        Jp = Jp';
        if V(i,j)~= 0
            x = [x;traj{j}(i,:)];
            [i1,i2] = meshgrid(2*j-1:2*j,7*i-6:7*i);
            J1 = sparse(i1(:),i2(:),Jc(:),2*Cam*n,7*Cam+3*n);
            [i1,i2] = meshgrid(2*j-1:2*j,7*Cam+3*j-2:7*Cam+3*j);
            J2 = sparse(i1(:),i2(:),Jp(:),2*Cam*n,7*Cam+3*n);
            J = J+J1+J2;
        else
            x = [x;zeros(1,2)];
            x_repro(1:2,j) = zeros(2,1);
        end
    end
    error_u = x(:,1)'-x_repro(1,:)./x_repro(3,:);
    error_v = x(:,2)'-x_repro(2,:)./x_repro(3,:);
    error = [error_u;error_v];
    e = [e;error(:)];
    
    
%     if nargout>1
%         for j = 1:n
%             Jc = [PDfq(K,R,C,X0(j,1:3)),PDfC(K,R,C,X0(j,1:3))];
%             Jc = Jc';
%             Jp = PDfX(K,R,C,X0(j,1:3));
%             Jp = Jp';
%             if V(i,j) ~= 0
%                 [i1,i2] = meshgrid(2*j-1:2*j,7*i-6:7*i);
%                 J1 = sparse(i1(:),i2(:),Jc(:),2*Cam*n,7*Cam+3*n);
%                 [i1,i2] = meshgrid(2*j-1:2*j,7*Cam+3*j-2:7*Cam+3*j);
%                 J2 = sparse(i1(:),i2(:),Jp(:),2*Cam*n,7*Cam+3*n);
%                 J = J+J1+J2;
%             end
%         end
%     end      
end
end