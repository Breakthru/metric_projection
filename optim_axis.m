function [out]=optim_axis(W,epsilon,n_iter,Minit,ind1,ind2,name)
% Articulated motion recovery using Metric Projection
%
% Ref: "Factorization for Non-Rigid and Articulated Structure using Metric Projections"
% Marco Paladini, Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% CVPR 2009, June 20-25 Miami, Florida
%
% Authors: Marco Paladini (paladini@dcs.qmul.ac.uk),  Alessio Del Bue, Marko Stošić, Marija Dodig, João Xavier, Lourdes Agapito
% Last Modified: 18/08/2009
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
% Starting of the powerfactorization cycle
while err_1(k)>epsilon && k<n_iter
  %--- the projection step
    tic();
    

    Motion=zeros(F2,5);  %this will satisfy motion constraints
    % now project the motion estimate
    for f=[1:2:F2]
        %% calculate u minimizing cost function       
        min_cf=inf;  % minimum of Cost Function
        x = Rprev(f:f+1,1);
        Y = Rprev(f:f+1,2:3);
        Z = Rprev(f:f+1,4:5);
        [Vy,Sy,Wy]=svd(Y);
        [Vz,Sz,Wz]=svd(Z);
        s1ps2 = Sy(1,1)^2 +Sy(2,2)^2;
        s1ms2 = Sy(1,1)^2 - Sy(2,2)^2;
        s1s2 = 2*Sy(1,1)*Sy(2,2);
        e1pe2 = Sz(1,1)^2 +Sz(2,2)^2;
        e1me2 = Sz(1,1)^2 - Sz(2,2)^2;
        e1e2 = 2*Sz(1,1)*Sz(2,2);
        % Compute cost function for every point in the grid
        bound=[-pi,pi,1e-5,1];
        %numpoints=0;
        %figure(); hold on
        for precision=[0.2,0.05,0.005,0.005,0.0005,0.0001,0.00001,5e-6]
        for R=(bound(3):precision:bound(4))
                Rsq = R^2;
              constRY = s1ps2 - Rsq*Sy(2,2)^2 + s1s2*sqrt(1-Rsq);
              constRZ = e1pe2 - Rsq*Sz(2,2)^2 + e1e2*sqrt(1-Rsq);
              % equally spaced points on a cirlce have radius:
              % -pi:precision/r:pi
              for th =(bound(1):precision/R:bound(2))
                v = [cos(th);sin(th)];
                u = R*v;
                %plot(u(1),u(2),'r.')
                cf1= - u(1)^2 - u(2)^2 - 2*u'*x;
                cf2 = - 2 * sqrt(constRY - Rsq*s1ms2*(v'*Vy(:,1))^2);
                cf3 = - 2 * sqrt(constRZ - Rsq*e1me2*(v'*Vz(:,1))^2);
                cf = cf1 + cf2 + cf3;
                %numpoints = numpoints+1;
                if cf < min_cf
                min_u = u;
                min_th=th;
                min_R=R;
                min_cf = cf;
                end
              end
           
        end
        bound=[min_th-precision/min_R,min_th+precision/R,min_R-precision,min_R+precision];
        if bound(4)>1 bound(4)=1; end
        %hold on
        %min_u
%         plot(min_u(1),min_u(2),'bo','LineWidth',4)
%         line([bound(3)*cos(bound(1)),bound(3)*cos(bound(2))],[bound(3)*sin(bound(1)),bound(3)*sin(bound(2))],'LineWidth',4)
%         line([bound(4)*cos(bound(1)),bound(4)*cos(bound(2))],[bound(4)*sin(bound(1)),bound(4)*sin(bound(2))],'LineWidth',4)
%         line([bound(3)*cos(bound(1)),bound(4)*cos(bound(1))],[bound(3)*sin(bound(1)),bound(4)*sin(bound(1))],'LineWidth',4)
%         line([bound(3)*cos(bound(2)),bound(4)*cos(bound(2))],[bound(3)*sin(bound(2)),bound(4)*sin(bound(2))],'LineWidth',4)
        end
        %keyboard

        %% Minimum of cost found: done calculating u
        A = get_a(min_u,Y);
        B = get_a(min_u,Z);
        % end of the projection
        Motion(f:f+1,:)=[min_u A B];
    end
timing=[timing toc()];

   %% powerfactorization step
%         Shape = pinv(Motion)*W;
%         R = W*pinv(Shape);
    %easy mod to powerfact step
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
out.timing  = timing ;    
end

function A=get_a(u,Y)
% returns A such that ||A-Y|| is minimized and u|A is stiefel
% find optimal A
P = sqrtm(eye(2)-u*u');
if ~ isreal(P) disp(u);error('P has imaginary');end
[U,D,V] = svd(P*Y);
Q = U*V';
A = P*Q;
end

