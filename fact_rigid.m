% Created by Alessio Del Bue and Lourdes Agapito
% Last Modified: 18/08/2009
% License: GPLv2
%
% Function: [S,M,T,L,scale] = fact_nonrigid(W, align)
%
% Factorization algorithm for rigid 3D shape recovering
%
%
% Factorization accept as input the measurement matrix W containing the  P
% pixels position in F frames. 
%
% The function calls:   fact_metric_constraint
%                       fact_align
%
% Input:
%
% W: measurement matrix of size 2F by P where F is # of frames and P # of
% feature points
% align: Reference frame for camera matrix rotation alignment.  
% 
% Output:
%
% S: shape matrix of dim 3 by P. Contains the 3D points
% M: motion matrix, it encodes the rotation parameters size 2F by 3
% T: translation vector, results from the centroid calculation. Size 2*F by
% 1
% scale: scaling value for all the coordinates

function [S,M,T,scale] = fact_rigid(W, align)

r = 3;
s = size(W);
T = W*ones(s(2),1)/s(2); % mean value for every x, y
 
W = W - T*ones(1,s(2)); % registered measurement matrix

Wx = W([1:2:s(1)],:);
maxx = max(max(abs(Wx)));
c = [0:2:s(1)]; c(1)=[];
Wy = W(c,:);
maxy = max(max(abs(Wy)));
scale = (maxx+maxy);

W=W/scale;


[O1,lambda,O2] = svd(W,0); % SVD factorize the registered measurement matrix W

singvalue = diag(lambda)';
singvalue = singvalue(1:r);

lambda1 = lambda(1:r,1:r);	% W has rank r = 3*K, we take the greatest singular values
rlambda1 = sqrt(lambda1);  

Mhat = O1(:,1:r) * rlambda1;  % Matrix contains rotation matrix and the configuration weights l
Shat = rlambda1 * O2(:,1:r)';   % Shape matrix contains the K shapes    

[M, S, Q] = fact_rigid_constraint(Mhat,Shat);
% 
% [M S]=fact_align(M,S,align);    % Align the motion matrix to the same reference frame

