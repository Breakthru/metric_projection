function [M,L,Rstief]=init_nonrigid(W,K)
% Construct initial motion matrix for the metric projection
% algorithm
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% W: Measurement matrix of size 2F by P where F is # of frames and P # of feature points 
% K: Number of basis shapes (rank of decomposition will be 3*K)
%
% Output:
%
% M: Motion matrix of size 2*F by 3K
% L: Deformation weights matrix, size F by K
% Rstief: Camera matrices stacked vertically, size 2*F by 3. Those are approximately Stiefel

M=[]; S=[]; L=[]; Rstief=[]; %prepare output

Kiter=1;

f= size(W,1)/2;

%start with rigid fact 
[S,Rstief]=fact_rigid(W,1);
%[Rstief,S]=rigidfact(W,1e-5,200);
% perform a rank3 svd
L=[];
L(:,1)=ones(f,1);
    W_r = W - Rstief*S; % residual 
    
for n=2:K


    [U,D,V]=svd(W_r);
    U=U(:,1:3);
    D=D(1:3,1:3);
    V=V(:,1:3)';

    M_r=U*sqrt(D);
    S_r=sqrt(D)*V;




    for i = 0:f-1

        Mf(1,:) = M_r(2*i+1,:);   % extract rowwise a 2 by 3*k submatrix from Mhat such that Mf = [l1*R1...lk*R1]
        Mf(2,:) = M_r(2*i+2,:);   % vectorize Mf such that [l1*r1,l1*r4,l1*r2,l1*r5...lk*r1,lk*r4,lk*r2,lk*r5,lk*r3,lk*r6]' 

            vec6Mf(1:3,1) = Mf(1,:)';
            vec6Mf(4:6,1) = Mf(2,:)';        

        L_ij = ([Rstief(1,:),Rstief(2,:)]*vec6Mf)'/2;
        L(i+1,n) = L_ij; % plug the F configuration weights in a sinlgle F by k matrix
    end
    
    W_r = W_r - M_r*S_r;
end

M=makeMfromRl(L,Rstief);

end
