function [M,S,t1,t2,scale,reprojection_error] = fact_OMP_articulated_FS(W,epsilon,n_iter,ind1,ind2) 
% Articulated motion recovery, hinge joint
% Orthographic camera, Metric Projection algorithm
% FS: Full Search, projection is done via brute force optimisation
% see paper for details
%
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
% W: Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% epsilon: Stopping parameter, (e.g. 1e-6). Iterations will end if ||R-Rprev|| is smaller than epsilon (R being the new estimate of motion matrix)
% n_iter: Maximum number of iterations (e.g. 10)
% ind1: Array of the indexes for the points in W belonging to the first object
% ind2: Array of the indexes for the points in W belonging to the second object
%
% Output:
%
%  M: Motion matrix of size 2*F by 5
%  S: Rank 5 shape matrix for the two objects, size 5 by P, where P is the total # of points
%  t1: 2D Centroid coordinates for first object
%  t2: 2D Centroid coordinates for second object
%  scale: Scale factor used to normalise input for numerical stability
%  reprojection_error: rms pixel error


% use Tresadern-Reid linear approximation as starting point
Minit = tr_axis(W,ind1,ind2);
% run the Metric Projection algorithm
out=optim_axis(W,epsilon,n_iter,Minit,ind1,ind2,name,gt3D);
M = out.Motion;
S = out.Shape;
t1 = out.t1;
t2 = out.t2;
scale = out.scale;
[ourS3D, axis3D,Wrepr] = makeS3dAxis(M,scale*S,t1,t2,ind1,ind2)
reprojection_error=norm([W(:,ind1) W(:,ind2)] - Wrepr,'fro')/norm(W,'fro')