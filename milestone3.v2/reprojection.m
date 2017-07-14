function Gerror = reprojection(Cset,Rset,X,K,x,y,V,idx_frame)

I = length(Cset);
J = size(X,1);

% Load images
for i = 1 : I
    k = idx_frame(i); % reproject in sequence
    str = sprintf('Milestone3_data/SfMProjectData_1/image%07d.bmp', k);
    img{i} = imread(str);
end

% J*I
x_real = x.*V;
y_real = y.*V;
x_rep = zeros(J,I);
y_rep = zeros(J,I);
num = 0;

for i = 1:I
    P = K*Rset{i}*[eye(3),-Cset{i}];
    for j = 1:J
        x_rep(j,i) = V(j,i)*(P(1,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]');
        y_rep(j,i) = V(j,i)*(P(2,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]');
        if V(j,i) ==1
            num = num+1;
        end
    end
end
error = [x_real-x_rep;y_real-y_rep];
Gerror = sum((reshape(error,[],1)).^2)/num;

for i = 1:I
    figure;
    imshow(img{i}); % reproject in sequence
    hold on
    for j = 1:J
        if V(j,i) == 1
        plot(x_real(j,i), y_real(j,i), 'gx','MarkerSize',5);
        plot(x_rep(j,i), y_rep(j,i), 'rx','MarkerSize',5);
        end
    end
    hold off
end

end