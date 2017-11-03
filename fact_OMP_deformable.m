function [M,S,Rstiefel,l,t,scale,reprojection_error] = fact_OMP_deformable(W,rank,epsilon,n_iter)
% Factorization for Structure From Motion
% Orthographic camera, Metric Projection Algorithm for deformable objects
% Ref: "Factorization for Non-Rigid and Articulated Structure using Metric Projections"
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% CVPR 2009, June 20-25 Miami, Florida
%
% Authors: Marco Paladini (paladini@dcs.qmul.ac.uk),  Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% Wn: Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% K: Number of basis shapes (rank of decomposition will be 3*K)
% epsilon: Stopping parameter, iterations will end if ||R-Rprev|| is smaller than epsilon (R being the new estimate of motion matrix)
% n_iter: Maximum number of iterations
%
% Output:
%
% M: Motion matrix of size 2*F by 3K
% Shape: Basis Shapes stacked vertically, size 3K by P where P is the number of points
% Rstiefel: Camera matrices stacked vertically, size 2*F by 3. Those are exactly Stiefel
% l: Deformation weights matrix, size F by K
% t: Recovered translation vector
% scale: Assigned scale to input data (input within -1:1 for numerical stability)
% reprojection_error: rms pixel error

if mod(rank,3)~=0
    error('basis shape model requires a rank multiple of 3')
end

K = rank/3;
[M,S,Rstiefel,l,err_1,timing,scale,t,reproj]=optim_nonrigid(W,K,epsilon,n_iter);
reprojection_error = norm(W - (M*S*scale + repmat(t',1,size(S,2))),'fro')/norm(W,'fro')
end
