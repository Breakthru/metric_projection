function [Z,T,err1,Shape,Motion,l,Rstief,timing,scale] = md_optim(Wo,D,K,tresh1,tresh2,maxIter1,maxIter2,name)
% Factorization for deformable shapes using metric projections with missing data.
% Extends Marques-Costeira [2] algorithm to deformbale structures
%
% Ref: "Factorization for Non-Rigid and Articulated Structure using Metric Projections" [1]
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% CVPR 2009, June 20-25 Miami, Florida
%
% Authors: Marco Paladini (paladini@dcs.qmul.ac.uk),  Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% Wo: Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% D: Visibility matrix, size 2F by P, 1 if point is visible 0 is missing
% K: Number of basis shapes (rank of decomposition will be 3*K)
% tresh1: Stopping parameter for missing data iterations
% tresh2: Stopping parameter for factorization. Iterations will end if ||R-Rprev|| is smaller than this (R being the new estimate of motion matrix)
% maxIter1: Maximum number of iterations for outer loop (see algorithm 2 in [1])
% maxIter2: Maximum number of iterations for inner loop (see algorithm 1 in [1])
% name: (for debugging) if specified, save all intermediate results in mat-files with prefix 'name'
%
% Output:
%
% M: Motion matrix of size 2*F by 3K
% Shape: Basis Shapes stacked vertically, size 3K by P where P is the number of points
% Rstief: Camera matrices stacked vertically, size 2*F by 3. Those are exactly Stiefel
% l: Deformation weights matrix, size F by K
% err_1: Vector of stopping condition values
% timing: Vector of time elapsed for each iteration (in seconds)
% scale: Assigned scale to input data (input within -1:1 for numerical stability)
% t: Recovered translation vector
% reproj: Calculated reprojection error for all iterations
%
% Ref [2]: "Estimating 3D shape from degenerate sequences with missing data"
% Manuel Marques and Joao Costeira
% Journal of Computer Vision and Image Understanding
% Volume 113 ,  Issue 2  (February 2009)


numPoints = size(Wo,2);

% scale=max(max(abs(Wo)));
% Wo=Wo/scale;
WoD=(Wo.*D);
tic();
    Wo(D==0)=nan;
    [Motion, Shape, Winit] = missingdata_rigid(Wo, D, 1e-7,100, 200);
    Wo(D==0)=0;
toc()
disp('done initialisation')
Zprev = Winit .* not(D) + Wo .* D;
%[u,m]=register(Z0.*D);
%%% put the mean values in the missing slots
%%% Z=Z0.*D + not(D).*(m'*ones(1,numPoints));
%Z = (Z0+randn(size(Z0)) ).*not(D) + Z0.*D;


k=1;
err1=[inf];
%normfact=max(max(abs(Zprev)));
%Zprev = Zprev/scale;

%tic();
timing=[];
while err1(end) > tresh1 && (k < maxIter1)
  [Z,T]=register(Zprev); % centroid estimation
  %step 3, estimate Motion and Shape
  
  [Motion,Shape,Rstief,l,err_1,time1,scale]=optim_nonrigid(Z,K,tresh2,maxIter2);
  % step 4, Update data matrix
  Z = ( Motion*(scale*Shape) + T'*ones(1,numPoints) ).* not(D)  + WoD;
  % reproject missing points + known points
  err1(k)= norm(Zprev - Z,'fro')/numel(Z);
  k = k+1;
  Zprev=Z;
  timing=[timing time1];
  if exist('name','var')
	save([name sprintf('_iter%02d.mat',k)],'Motion','Shape','Rstief','l','time1','Z','Wo','D','K','Winit','T');
  end

end



end
