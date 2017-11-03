function [M,S,Zout]=missingdata_rigid(Z,D,epsilon,n_iter,n_iter2)
% Rigid Factorization with missing data
% Author: Marco Paladini
% Last Modified: 12/10/2009
% License: GPLv2
%
% Ref: "Estimating 3D shape from degenerate sequences with missing data"
% Manuel Marques and Joao Costeira
% Journal of Computer Vision and Image Understanding
% Volume 113 ,  Issue 2  (February 2009)
% 
% Special thanks to Manuel Marques for finding bugs and correcting
% my implementation of his algorithm
%
% Input:
% Z:  Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% D: Visibility matrix, size 2F by P, 1 if point is visible 0 is missing
% tresh1: Stopping parameter for iterations
% n_iter: Maximum number of iterations for missing data estimation
% n_iter2: Maximum number of iterations for rigid factorization
%
% Output:
%
% M: Motion matrix, rank 3
% S: Shape matrix, rank 3
% Zout: Measurement matrix with the estimated missing data.

[F2,P]=size(Z);
Dhat=ones(F2,P)-D;
Zk=Z;
k=1;
Z(D == 0) = 0;
Zk = (sum(Z,2)./sum(D,2)*ones(1,size(Z,2))).*not(D) + Z.*D;
ref_t = find(sum(D,1) == F2, 1);
num_missing = sum(sum(not(D)));

if ~isempty(ref_t)
    t = Z(:,ref_t(1))';
    Zc = Zk - repmat(t',1,P);
else
    [Zc,t]=register(Zk);
end
[U,sv,V]=svd(Zk);
sv=sv(1:3,1:3);
U=U(:,1:3);
R=U*sqrtm(sv);
% if there is no missing data return a metric constrained factorization
if num_missing==0
  warning('No missing data, using full data for rigid factorization...');
  [M,S,R]=rigidfact(Zc,epsilon,n_iter2,R);
  Zout=Zc+ t'*ones(1,P); % no missing data to re-estimate
  return
end

thenorm=[inf];
while thenorm(k)>epsilon && k<=n_iter
    Zkp=Zc;
    [Mk,Sk,R]=rigidfact(Zc,epsilon,n_iter2,R);
    % update data matrix
    Zk=Mk*Sk.*Dhat + Zc.*D + t'*ones(1,P) ;
    k=k+1;
    thenorm=[thenorm sum(sum(abs(Zk-Zkp-repmat(t',1,P))))/num_missing]
    if ~isempty(ref_t)
        t = Z(:,ref_t(1))';
        Zc = Zk - repmat(t',1,P);
    else
        [Zc,t]=register(Zk);
    end
%    disp([num2str(k) ' steps performed'])
end
M=Mk;
S=Sk;
Zout=Zk;
end

function [M,S,R]=rigidfact(Zc,epsilon,n_iter,R)
% Rigid Factorization
% Marco Paladini, 15/03/2008
% reference: Optimal Shape from motion
% estimation with missing and degenerate data
% Manuel Marques and Joao Costeira

% 1. Factorize Z_c using any factorization
% 2. project R into the manifold of motion matrices
% 3. recalculate shape with S^k = {M^k}^\daga Z_c
% 4. recalculate R as Z_c {S_k}^\daga
% 5. if ||M_k - M_{k-1} || < \epsilon finish
% 6. M=M_k S=S_k

if nargin < 2
    epsilon=0.05
end

F=size(Zc,1);
k=0;
thenorm=epsilon+1;
while thenorm>epsilon && k<n_iter
    Rp=R;
    M=[];
%     for f=1:2:F
%         Rf=R(f:f+1,:);
%         [U,S,V]=svd(Rf);
%         U(:,3)=[0;0];
%         U(3,:)=[0 0 0];
%         T=U*V;
%         M=[M;T(1:2,:)];
%     end
    for f=1:2:F
       Rf=R(f:f+1,:);
       [U,S,V]=svd(Rf','econ');
       T=mean(diag(S))*U*V';
       M=[M;T'];
   end 
    S=pinv(M)*Zc;
    R=Zc*pinv(S);
    k=k+1;
    thenorm=sum(sum(abs(R-Rp)))/(F/2);
    %disp([num2str(k) ' steps performed'])
end

end