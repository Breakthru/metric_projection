function [Minit,STR,t1,t2,normfact]=tr_axis(W,ind1,ind2)
% Articulated Hinge joint recovery by Factorization
% This is an implementation of Tresadern-Reid method for hinge
% joint recovery, see reference for further details
%
% Author: Marco Paladini, João Fayad, Alessio Del Bue, Lourdes Agapito
% Last Modified: 18/08/2009
%
% Ref: "Articulated Structure from Motion by Factorization"
% Tresadern, P. and Reid, I.
% CVPR '05: Proceedings of the 2005 IEEE Computer Society Conference on Computer Vision and Pattern Recognition
% 
% Input:
% W: 2d tracks of feature points, each frame stacked vertically, each
% column represents a point
% ind1: Array of the indexes for the points in W belonging to the first object
% ind2: Array of the indexes for the points in W belonging to the second
% object
%
% Output:
% Minit: Resulting Motion matrix
% STR: Combined Shape matrix of the two objects, rank 5
% t1: Centroid coordinates for first object
% t2: Centroid coordinates for second object
% normfact: Scaling factor applied to input data for numerical stability


        W1=W(:,ind1);
        W2=W(:,ind2);
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
       	normfact=max(max(abs(W)));
        W=W/normfact;

        [U,D,V]=svd(W);
        % Impose block structure

        D = D(1:5,1:5);
        U = U(:,1:5);
        V = V(:,1:5)';

        W = U*D*V;
        ind1=[1:size(W1,2)];
        ind2=[size(W1,2)+1:size(W,2)];
        [Minit,S_]=block_axis(U*D,V,ind1,ind2,1);

        % approximate metric upgrade:
        Q1=make_q(Minit(:,1:3));
        Q2=make_q(Minit(:,[1 4 5]));
        Q2=Q2* (Q1(1,1)/Q2(1,1));
        Q=[[Q1;zeros(2,3)] [Q2(1,2:3);zeros(2,2);Q2(2:3,2:3)]];
        %Q=Q/Q(5,5);
        Minit=Minit*Q;
        STR=inv(Q)*S_;
        
end

function Q= make_q(M)
% Find linear transformation Q
% such that M*Q has the rows pairwise orthogonal


F = size(M,1)/2;

A = zeros(F * 3,6);
b = zeros(1,F * 3);

for n=[1:2:2*F]
    vec_i(n,:)=M(n,:);
    vec_j(n,:)=M(n+1,:);
end

% compute the 3Fx6 matrix that represent the system of metric constraits
% as Aq=y to find the elements q of the matrix Q
for k=[1:F]
        A(k,:)=[vec_i(k,1)^2 2*vec_i(k,1)*vec_i(k,2) 2*vec_i(k,1)*vec_i(k,3) vec_i(k,2)^2 2*vec_i(k,2)*vec_i(k,3) vec_i(k,3)^2];
        A(k+F,:)=[vec_j(k,1)^2 2*vec_j(k,1)*vec_j(k,2) 2*vec_j(k,1)*vec_j(k,3) vec_j(k,2)^2 2*vec_j(k,2)*vec_j(k,3) vec_j(k,3)^2];
        A(k+2*F,:)=[vec_i(k,1)*vec_j(k,1) vec_i(k,2)*vec_j(k,1)+vec_i(k,1)*vec_j(k,2) vec_i(k,3)*vec_j(k,1)+vec_i(k,1)*vec_j(k,3) vec_i(k,2)*vec_j(k,2) vec_i(k,3)*vec_j(k,2)+vec_i(k,2)*vec_j(k,3) vec_i(k,3)*vec_j(k,3)];
end

% finished building A matrix (3Fx9)
% build b
for k=[1:2*F]    b(k)=1; end
for k=[2*F+1:3*F] b(k)=0; end
% compute Q matrix
q=pinv(A)*b';
Q=[q(1) q(2) q(3);q(2) q(4) q(5);q(3) q(5) q(6)];
Q=inv(chol(inv(Q)));


end

function v = zt(a, b)

v = [    a(1)*b(1), a(1)*b(2)+a(2)*b(1), a(1)*b(3)+a(3)*b(1), ...
   a(2)*b(2), a(2)*b(3)+a(3)*b(2), a(3)*b(3) ];

end

function [M,S]=block_axis(U,V,ind1,ind2,axis)
% Applying the Nl(.) operator described in the article. The svd computation
% of the null space is used to avoid the empty null space caused by
% numerical errors.

V1 = V(:,ind1);
[A,B,C] = svd(V1);
NL1 = A(:,end-1:end);

V2 = V(:,ind2);
[A,B,C] = svd(V2);
NL2 = A(:,end-1:end);

% Defining the Au matrix described in the article
if axis==1
    Au = [ 1 0 0 0 0;NL2' ; NL1']; %x axis
else if axis==2
    Au = [ NL2(:,1)' ;0 1 0 0 0; NL2(:,2)'; NL1']; %y axis
    else if axis==3
Au = [ NL2' ;0 0 1 0 0; NL1']; %z axis
        end
    end
end



% Recovering the motion M and object S matrices:
S = Au * V;
M =  U * inv(Au);
end
