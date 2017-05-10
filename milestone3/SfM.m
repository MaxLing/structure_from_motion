%% load data
 sequence = [1 2 3 4 5 6];
% sequence = [1 4 2 3 5 6];
% sequence = [3 4 1 2 5 6];

load DataRANSAC0.5.mat

%%  Triangulation
initialimg1=sequence(1);  %23 34 45
initialimg2=sequence(2);

idx = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-initialimg1))+initialimg2-initialimg1)==1);
x1 = [MatchX(idx,initialimg1) MatchY(idx,initialimg1)];
x2 = [MatchX(idx,initialimg2) MatchY(idx,initialimg2)];

F = EstimateFundamentalMatrix(x1, x2); % Estimate F and E.
E = EssentialMatrixFromFundamentalMatrix(F, K);   
[Cs, Rs] = ExtractCameraPose(E); % Extract C and R.
% Triangulate Points for each C and R.
for i = 1:4
    Xs{i} = LinearTriangulation(K, zeros(3,1), eye(3), Cs{i}, Rs{i}, x1, x2);
end
% Disambiguate the actual pose.
[C, R, X] = DisambiguateCameraPose(Cs, Rs, Xs);
% Gerror1 = reprojection({zeros(3,1),C},{eye(3),R},X,K,MatchX(idx,idx_frame),MatchY(idx,idx_frame),ones(length(idx),2),idx_frame);
%PlotCamerasAndPoints( {zeros(3,1),C}, {eye(3),R}, X, 1 );
% remove weird points
idx_error = find(X(:,3) < 0 | sqrt(sum(X.^2,2)) > 80);
X(idx_error,:) = [];
idx(idx_error,:) = [];
x1(idx_error,:) = [];
x2(idx_error,:) = [];
% Gerror2 = reprojection({zeros(3,1),C},{eye(3),R},X,K,MatchX(idx,idx_frame),MatchY(idx,idx_frame),ones(length(idx),2),idx_frame);
%PlotCamerasAndPoints( {zeros(3,1),C}, {eye(3),R}, X, 1 );

% nonlinear
X = NonlinearTriangulation(K,zeros(3,1),eye(3),C,R,x1,x2,X);
% Gerror3 = reprojection({zeros(3,1),C},{eye(3),R},X,K,MatchX(idx,idx_frame),MatchY(idx,idx_frame),ones(length(idx),2),idx_frame);
%PlotCamerasAndPoints( {zeros(3,1),C}, {eye(3),R}, X, 1 );
% remove weird points
idx_error = find(X(:,3) < 0 | sqrt(sum(X.^2,2)) > 80);
X(idx_error,:) = [];
idx(idx_error,:) = [];
% Gerror4 = reprojection({zeros(3,1),C},{eye(3),R},X,K,MatchX(idx,idx_frame),MatchY(idx,idx_frame),ones(length(idx),2),idx_frame);
%PlotCamerasAndPoints( {zeros(3,1),C}, {eye(3),R}, X, 1 );

idx_frame = [initialimg1,initialimg2]; % index for frame

Cset{1} = zeros(3,1);
Rset{1} = eye(3);
Cset{2} = C;
Rset{2} = R;

X3D(idx,:) = X; % idx here is visible 3D point
V(idx, initialimg1) = 1;
V(idx, initialimg2) = 1;

%% PnP
for s = 3:img_num % register new cameras from 3rd
    j = sequence(s);
        idx_corr = [];
        for k = 1:length(idx_frame)
            if idx_frame(k) < j
            idx_old = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-idx_frame(k)))+j-idx_frame(k))==1);
            else
            idx_old = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-j))+idx_frame(k)-j)==1);
            end
            idx_corr = union(idx_corr,idx_old);
        end
        idx3 = intersect(idx,idx_corr);
 
% %       %visualization        
%         V(idx3, j) = 1;
%         X = X3D(idx3,:);
%         x = [MatchX(idx3,j),MatchY(idx3,j)];
%         [Cnew,Rnew] = LinearPnP(X,x,K);
%         Gerror1 = reprojection({Cnew},{Rnew},X,K,MatchX(idx3,j),MatchY(idx3,j),V(idx3, j),j);
%         PlotCamerasAndPoints( {zeros(3,1),C,Cnew}, {eye(3),R,Rnew}, X, 1 );
%         
%         [X,x,inliners_id] = PnPRANSAC(X,x,K);
%         idx3 = idx3(inliners_id);
%         [Cnew,Rnew] = LinearPnP(X,x,K);
%         Gerror2 = reprojection({Cnew},{Rnew},X,K,MatchX(idx3,j),MatchY(idx3,j),V(idx3, j),j);
%         PlotCamerasAndPoints( {zeros(3,1),C,Cnew}, {eye(3),R,Rnew}, X, 1 );
%         
%         [Cnew,Rnew] = NonlinearPnP(X,x,K,Cnew,Rnew);
%         Gerror3 = reprojection({Cnew},{Rnew},X,K,MatchX(idx3,j),MatchY(idx3,j),V(idx3, j),j);
%         PlotCamerasAndPoints( {zeros(3,1),C,Cnew}, {eye(3),R,Rnew}, X, 1 );
        
        if isempty(idx3)
            continue;
        elseif length(idx3) < 6
            X = X3D(idx3,:);
            x = [MatchX(idx3,j),MatchY(idx3,j)];
            [Cnew,Rnew] = LinearPnP(X,x,K);
            [Cnew,Rnew] = NonlinearPnP(X,x,K,Cnew,Rnew);
%         if length(idx3) < 6
%             continue;
        else
        X = X3D(idx3,:);
        x = [MatchX(idx3,j),MatchY(idx3,j)];
        
        [X,x,inliners_id] = PnPRANSAC(X,x,K);
        idx3 = idx3(inliners_id);
        [Cnew,Rnew] = LinearPnP(X,x,K);
        [Cnew,Rnew] = NonlinearPnP(X,x,K,Cnew,Rnew);  % register ith image
        end
        Cset{end+1} = Cnew;
        Rset{end+1} = Rnew;
        idx_frame(end+1) = j;
        V(idx3, j) = 1;
        
%% add 3d points in current frames
        for i = 1:length(idx_frame)-1
            idx1 = setxor(idx,(1:size(Match,1))'); % new points
            if idx_frame(i) < j
            idx2 = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-idx_frame(i)))+j-idx_frame(i))==1);
            else
            idx2 = find(Match(:,sum(1:(img_num-1))-sum(1:(img_num-j))+idx_frame(i)-j)==1);
            end
%             idx2 = find(Match(:,j)==1 & Match(:,idx_frame(i))==1); % visible in previous and new image (and not outliners)
            idx4 = intersect(idx1,idx2); % idx4 here is new 3D points
            x1 = [MatchX(idx4,idx_frame(i)), MatchY(idx4,idx_frame(i))];
            x2 = [MatchX(idx4,j), MatchY(idx4,j)];
            
            Xnew = LinearTriangulation(K, Cset{i}, Rset{i}, Cnew,Rnew, x1, x2);
            % remove weird points
            idx_error = find(Xnew(:,3) < 0 | sqrt(sum(Xnew.^2,2)) > 80);
            Xnew(idx_error,:) = [];
            idx4(idx_error,:) = [];
            x1(idx_error,:) = [];
            x2(idx_error,:) = [];
            
            Xnew = NonlinearTriangulation(K, Cset{i}, Rset{i}, Cnew,Rnew, x1, x2,Xnew); % X0 is initial estimate from linear triangulation
            % remove weird points
            idx_error = find(Xnew(:,3) < 0 | sqrt(sum(Xnew.^2,2)) > 80);
            Xnew(idx_error,:) = [];
            idx4(idx_error,:) = [];
            
            % add new points, update X3D, idx, V
            idx = union(idx,idx4);
            X3D(idx4,:) = Xnew;
            V(idx4,j) = 1;
            V(idx4,idx_frame(i)) = 1;
        end
end

%% Bundle adjustment - uncomment to use
% first use idx to extract reconstructed points 
% [Cset, Rset, Xset] = BundleAdjustmentWithJacobian(Cset,Rset,X3D(idx,:),K, MatchX(idx,idx_frame),MatchY(idx,idx_frame),V(idx,idx_frame)); 

Xset = X3D(idx,:);

%% Visualization
PlotCamerasAndPoints( Cset, Rset, Xset, 1 );
Gerror = reprojection(Cset,Rset,Xset,K,MatchX(idx,idx_frame),MatchY(idx,idx_frame),V(idx,idx_frame),idx_frame);