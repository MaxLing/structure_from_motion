function gerror = TestReprojection(C,R,X,K,x,y,V,t)

J = size(X,1);
x = x(:,t);
y = y(:,t);
V = V(:,t);

% J*I
x_real = x.*V;
y_real = y.*V;
x_rep = zeros(J,1);
y_rep = zeros(J,1);
num = 0;

    P = K*R*[eye(3),-C];
    for j = 1:J
        x_rep(j) = V(j)*(P(1,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]');
        y_rep(j) = V(j)*(P(2,:)*[X(j,:) 1]')/(P(3,:)*[X(j,:) 1]');
        if V(j) ==1
            num = num+1;
        end
    end

error = [x_real-x_rep;y_real-y_rep];
gerror = sum((reshape(error,[],1)).^2)/num;

%     figure;
%     str = sprintf('Milestone3_data/SfMProjectData_1/image%07d.bmp', t);
%     img{t} = imread(str);
%     imshow(img{t}); % reproject in sequence
%     hold on
%     for j = 1:J
%         if V(j) == 1
%         plot(x_real(j), y_real(j), 'gx','MarkerSize',5);
%         plot(x_rep(j), y_rep(j), 'rx','MarkerSize',5);
%         end
%     end
%     hold off

end