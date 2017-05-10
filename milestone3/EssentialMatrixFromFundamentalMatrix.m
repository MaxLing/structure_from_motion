function [E] = EssentialMatrixFromFundamentalMatrix(F,K)

% Input: F K1 K2
% Output: E

% E = K2'*F*K1; 
E = K'*F*K; 

% 2nd cleanup
[U,D,V]=svd(E);
D = [1 0 0;0 1 0;0 0 0];
E = U*D*V';

end