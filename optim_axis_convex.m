function [out]=optim_axis_convex(W,epsilon,n_iter,Minit,ind1,ind2,name)
% Articulated motion recovery using Metric Projection
%
% Ref: "Metric Projections for Deformable and Articulated Structure-From-Motion"
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% Journal paper in preparation
%
% Authors: Marco Paladini (paladini@dcs.qmul.ac.uk),  Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% Last Modified: 6/10/2009
% License: GPLv2
%
% Reference for Articulated motion model: 
% "Articulated Structure from Motion by Factorization"
% Tresadern, P. and Reid, I.
% CVPR '05: Proceedings of the 2005 IEEE Computer Society Conference on Computer Vision and Pattern Recognition
%
% Input:
%
% W: Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% epsilon: Stopping parameter, iterations will end if ||R-Rprev|| is smaller than epsilon (R being the new estimate of motion matrix)
% n_iter: Maximum number of iterations
% Minit: Initial estimate of the motion matrix
% ind1: Array of the indexes for the points in W belonging to the first object
% ind2: Array of the indexes for the points in W belonging to the second object
% name: (for debugging) if specified, save all intermediate results in mat-files with prefix 'name'
%
% Output:
%
%  out.Motion: Motion matrix of size 2*F by 5
%  out.Shape: Rank 5 shape matrix for the two objects, size 5 by P, where P is the total # of points
%  out.t1: Centroid coordinates for first object
%  out.t2: Centroid coordinates for second object
%  out.scale: Scale factor used to normalise input for numerical stability
%  out.err_1: Vector of stopping condition values
%  out.timing: Vector of time elapsed for each iteration (in seconds)


W1 = W(:,ind1);
W2 = W(:,ind2);
[U,D,V]=svd(W1);
U=U(:,1:4);
D=D(1:4,1:4);
V=V(:,1:4)';
W1=U*D*V;
[U,D,V]=svd(W2);
U=U(:,1:4);
D=D(1:4,1:4);
V=V(:,1:4)';
W2=U*D*V;
[W1,t1]=register(W1);
[W2,t2]=register(W2);
W=[W1 W2];
scale=max(max(abs(W)));
W=W/scale;
[U,D,V]=svd(W);

D = D(1:5,1:5);
U = U(:,1:5);
V = V(:,1:5)';

W = U*D*V;
W1 = W(:,ind1);
W2 = W(:,ind2);
F2 = size(W,1);
k=1;
err_1=[inf];

Rprev=Minit;
d=[]; for f=1:2:F2 d=[d norm(Rprev(f:f+1,1))]; end
d = max(d);
Rprev = Rprev / d;
scale = d*scale;
    vals=[]; maxval=[];
timing=[];
%% Starting of the powerfactorization cycle
while err_1(k)>epsilon && k<n_iter
  %--- the projection step
    tic();
    Motion=zeros(F2,5);  %this will satisfy motion constraints
    num_frames = F2/2; % number of frames
    % now project the motion estimate
progress=0; % percentage
for i = 1:num_frames
    
    % print periodically the current problem number
    if rem(i,floor(num_frames/5)) == 0; progress=progress+20; fprintf('%d%%... ',progress); end;
    % retrieve problem instance from the stack
    x = Rprev((i-1)*2+1:i*2,1);
    Y = Rprev((i-1)*2+1:i*2,2:3);
    Z = Rprev((i-1)*2+1:i*2,4:5);
        [u,f]=relaxation_articulated(x,Y,Z);
        A = get_a(u,Y);
        B = get_a(u,Z);
        % end of the projection
        Motion((i-1)*2+1:i*2,:)=[u A B];
   
end;

        timing=[timing toc()];
   %% powerfactorization step
%         Shape = pinv(Motion)*W;
%         R = W*pinv(Shape);
%	  Pseudo-inversion of the two objects is done separately
    S1=pinv(Motion(:,[1 2 3]))*W(:,ind1);
    S2=pinv(Motion(:,[1 4 5]))*W(:,ind2);
    Shape=[S1(1,:) S2(1,:);
        S1(2,:) zeros(size(ind2));
        S1(3,:) zeros(size(ind2));
        zeros(size(ind1)) S2(2,:);
        zeros(size(ind1)) S2(3,:)];
    R = W*pinv(Shape);
    
    disp([num2str(k) ' iterations done.'])
    k=k+1;
    err_1(k)=norm(R-Rprev,'fro')/numel(R);
    Rprev=R;
    if exist('name','var')
        matname = [name '_' num2str(n_iter) 'iter_step_' num2str(k) '_axis.mat'];
        save(matname,'Motion','Shape','t1','t2','scale','W');
    end
end
out.Motion = Motion;
out.Shape = Shape;
out.t1 = t1;
out.t2 = t2;
out.scale = scale;
out.err_1 = err_1;
out.timing  = timing;
end

function A=get_a(u,Y)
% returns A such that ||A-Y|| is minimized and u|A is stiefel
% find optimal A
P = sqrtm(eye(2)-u*u');
if ~ isreal(P) disp(u);error('p has imaginary');end
[U,D,V] = svd(P*Y);
Q = U*V';
A = P*Q;
end


function [u_relaxation,objective_score] = relaxation_articulated(x,Y,Z)
cvx_quiet(true);

% compute constants
R1 = Y*Y';
n1 = trace(R1);
d1 = abs(det(Y));

R2 = Z*Z';
n2 = trace(R2);
d2 = abs(det(Z));

I2 = eye(2);

aux = pinv([ I2(:) R1(:) R2(:) ])*kron(x,x); 
alfa = aux(1);
beta = aux(2);
gamma = aux(3);

% solve convex relaxation
cvx_begin
    
    variable U(2,2) symmetric;
        
    maximize(trace(U)+2*sqrt(alfa*trace(U)+beta*trace(R1*U)+gamma*trace(R2*U))+2*sqrt(n1-trace(R1*U)+2*d1*sqrt(1-trace(U)))+2*sqrt(n2-trace(R2*U)+2*d2*sqrt(1-trace(U))));
        
    U == semidefinite(2);
    trace(U) <= 1;
        
cvx_end;
[Uaux,Saux,Vaux] = svd(U);
u_relaxation = Uaux(:,1)*sqrt(Saux(1,1));
if x'*u_relaxation < 0
    u_relaxation = -u_relaxation;
end;



% extract solution
U = u_relaxation*u_relaxation';
objective_score = trace(U)+2*sqrt(alfa*trace(U)+beta*trace(R1*U)+gamma*trace(R2*U))+2*sqrt(n1-trace(R1*U)+2*d1*sqrt(1-trace(U)))+2*sqrt(n2-trace(R2*U)+2*d2*sqrt(1-trace(U)));
end