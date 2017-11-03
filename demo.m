% Demo program for the metric_projection package,
% Ref: "Factorization for Non-Rigid and Articulated Structure using Metric Projections"
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% CVPR 2009, June 20-25 Miami, Florida
%
%  This program is free software; you can redistribute it and/or
%  modify it under the terms of the GNU General Public License
%  as published by the Free Software Foundation; version 2, June 1991 
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% If you use this code for a scientific publication, please reference the original paper:
% "Factorization for Non-Rigid and Articulated Structure using Metric Projections"
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% CVPR 2009, June 20-25 Miami, Florida
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% Output:
%
% load the synthetic sequence
load data.mat GT W2d
numframes = 74; % number of frames to work with
K = 5; % number of basis shapes, rank = 3*K

disp ( ['using ', num2str(numframes) ,' frames'])
W2d = W2d(1:2*numframes,:);
GT.shape3D = GT.shape3D(1:3*numframes,:);




% run the metric projection method
disp('running metric projection...')
[M,S,Rstiefel,l,t,scale,reprojection_error] = fact_OMP_deformable(W2d,3*K,1e-7,15);
% Basis shapes and coefficients make the reconstructed shape
rec_3D = makeS3fxp(l,S);
% show the results

[val,our3D,gt3D] = check_multiple_shapes(rec_3D,GT.shape3D);

disp(['3d error: ' num2str(mean(val)*100) ' % '])
visualise(W2d,gt3D,our3D)

