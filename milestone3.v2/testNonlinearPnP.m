[Cnew1,Rnew1] = NonlinearPnP(X,x,K,Cnew,Rnew);

Cset1{1} = zeros(3,1);
Rset1{1} = eye(3);
Cset1{2} = Cnew;
Rset1{2} = Rnew;
PlotCamerasAndPoints( Cset1, Rset1, X, 1 );

Cset2{1} = zeros(3,1);
Rset2{1} = eye(3);
Cset2{2} = Cnew1;
Rset2{2} = Rnew1;
PlotCamerasAndPoints( Cset2, Rset2, X, 1 );