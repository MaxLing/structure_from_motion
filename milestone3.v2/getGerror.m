function [gerror,error] = getGerror(Cset,Rset,X,K,x,y,V)

I = length(Cset);
J = size(X,1);
error = zeros(I*J*2,1);
num = 0;

for i = 1:I
    P = K*Rset{i}*[eye(3),-Cset{i}];
    for j = 1:J
        error(((i-1)*J+j-1)*2+1:((i-1)*J+j-1)*2+2) = [V(j,i)*(x(j,i) - (P(1,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]'));...
                                                        V(j,i)*(y(j,i) - (P(2,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]'))];
        if V(j,i) ==1
            num = num+1;
        end
    end
end

    gerror = sum(error.^2)/num;
end