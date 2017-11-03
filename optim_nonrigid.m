function [M,Shape,Rstief,l,err_1,timing,scale,t,reproj]=optim_nonrigid(Wn,K,epsilon,n_iter,name)
% Factorization for deformable shapes using metric projections
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



%disp(['using ' num2str(K) ' basis shapes'])
W2d = Wn;
[Wn,t]=register(Wn);
scale=max(max(abs(Wn)));
Wn=Wn/scale;

F2 = size(Wn,1);
k=1;
err_1=[inf];
reproj=[inf];
err3d=[inf];
Rprev=init_nonrigid(Wn,K);
timing=[];
% Starting of the powerfactorization cycle
best.k=1;
best.reproj=inf;
can_stop=0;
while err_1(k)>epsilon && ~can_stop
  %--- the projection step
    tic();
    Rstief=zeros(F2,3); l=zeros(F2/2,K);
    M=zeros(F2,3*K);
    [Rstief,l]=projectNonrigid(Rprev);
    M = makeMfromRl(l,Rstief);
    Shape = (inv(M'*M)*M')*Wn;
    R = Wn*(Shape'*inv(Shape*Shape'));
    k=k+1;
    err_1(k)=norm(R-Rprev,'fro')/numel(R);
    reproj(k)=norm(W2d - (M*Shape*scale + repmat(t',1,size(Shape,2))),'fro')/norm(W2d,'fro');
    % backtrace if a better reprojection error was found earlier
    if reproj(k)<best.reproj
        best.reproj=reproj(k);
        best.k=k;
        best.l=l;
        best.Rstief=Rstief;
        best.Shape=Shape;   
    end
    if k>8
        check = reproj(k-8:k-1) < reproj(k-7:k);
        if sum(check)==numel(check) || k>n_iter
        can_stop=1;
        else
        can_stop=0;
        end
    end
    if exist('GT','var')
        [ourS3D] = makeS3fxp(l,Shape);
        val= check_multiple_shapes(ourS3D,GT.shape3D);
        err3d(k)=mean(val);
    end
    disp (['Iter ' num2str(k)]);
    Rprev=R;
    timing(k)=toc();

%save results
    if exist('name','var') && size(name,2)>0
      matname = [name '_' num2str(n_iter) 'iter_step_' num2str(k) '_' num2str(K) 'basis.mat'];
      save(matname,'Rstief','l','Shape','K','err_1','timing','Wn','scale','t','reproj');
    end
end
disp(['best k was: ' num2str(best.k)])
Shape=best.Shape;
Rstief=best.Rstief;
l=best.l;
M=makeMfromRl(l,Rstief);
end