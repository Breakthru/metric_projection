% Created by Alessio Del Bue and Lourdes Agapito
% Last Modified: 18/08/2009
% License: GPLv2
%
% Function: [M, S] = fact_metric_constraint(Mhat,Shat)
%
% Forces the metric constraint k times for each columnwise triplets in the
% 2Fx3k motion matrix M, it applies also the inverse transformation on the
% structure matrix S to maintain the coherence between the data. Based on
% the orthonormal normalization showed on Matthew's papers... .It extends
% the BHB approach (no 2nd svd of rank-1)
%
% Input
%
% Mhat: 2Fx3k motion matrix resulting from the 1st svd
% Shat: 3kxP structure matrix resulting from the 1st svd
%
% Output
%
% M: 2Fx3k metric constrained motion matrix 
% S: 3kxP structure matrix with inverse transformation

function [M, S, Q] = fact_rigid_constraint(Mhat,Shat)

k = size(Mhat, 2)/3;
f = size(Mhat,1)/2; % prior rigid model
M = zeros(2*f, 3*k);
S = zeros(size(Shat));
M3 = zeros(2*f, 3);
Q = zeros(3*k,3*k);

for i=0:k-1
    M3 = Mhat(:,3*i +1: 3*i+3);   % extract the k-th 2*F by 3 matrix (columns-wise)    

% Build the Q transformation matrix for M3. Q forces the orthogonality and the equality
% of the norm for every k 2*F by 3 matrix

	for n = 0:f-1,
      A(n+1,:) = fact_leq(M3(2*n+1,:), M3(2*n+1,:)); 
      A(n+1+f,:) = fact_leq(M3(2*n+2,:), M3(2*n+2,:)); % equal norm condition
      A(n+1+2*f,:) = fact_leq(M3(2*n+1,:), M3(2*n+2,:));   % orthogonality condition
%       A(n+1+2*f,:) = fact_leq(M3(2*n+2,:), M3(2*n+1,:));   % orthogonality condition
	end
% 	A(2*f+1,:) = fact_leq(M3(1,:), M3(1,:)); % fixed scale condition
 
% Solve for v, the 6 unknown values of the symmetric value
	b = [ones(2*f,1);zeros(f,1)];
	v = A \ b;    % Q*v = b
    
% C is a symmetric 3x3 matrix such that C = G * G'

	C(1,1) = v(1);                  
	C(1,2) = v(2);                  
	C(1,3) = v(3);                  
	C(2,2) = v(4);
	C(2,3) = v(5);
	C(3,3) = v(6);
	C(2,1) = C(1,2);
	C(3,1) = C(1,3);
	C(3,2) = C(2,3);
    
% 	[V D] = eig(C,'nobalance'); % not working, negatives eig
% 	G = V*sqrt(D);

% 	Qk = chol(C);
%   Q(3*i+1:3*i+3,3*i+1:3*i+3) = Qk;
    
%     
    [O1,lambda,O2] = svd(C);
    G = O1*sqrt(lambda);
	for m = 0:f-1
        num(2*m+1:2*m+2, :)  = M3(2*m+1:2*m+2, :)*G;
        den(2*m+1:2*m+2, :)  = pinv((M3(2*m+1:2*m+2, :)*G)');
	end
	Qk = G*sqrtm(num\den);
    Q(3*i+1:3*i+3,3*i+1:3*i+3) = Qk;
end

M = Mhat*Q;
S = inv(Q)*Shat;

end

% Created by Alessio Del Bue and Lourdes Agapito
% Last Modified: 18/08/2009
% License: GPLv2
%
% Function: v = fact_leq(r1, r2)
%
% Creates equations ( 1 by 6 vector) for the system computing the elements of
% the symmetric matrix
%
% Input:
%
% r1: 1st row of the rotation matrix
% r2: 2nd row of the rotation matrix
%
% Output:
%
% v = 6-vector whose elements compose the system values


function v = fact_leq(r1, r2)

v = [r1(1)*r2(1), r1(1)*r2(2)+r1(2)*r2(1), r1(1)*r2(3)+r1(3)*r2(1), ...
	r1(2)*r2(2), r1(2)*r2(3)+r1(3)*r2(2), r1(3)*r2(3) ];
end