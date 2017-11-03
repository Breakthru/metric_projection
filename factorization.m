
function [S,M,L] = factorization(Y,Sinitial,Minitial)

% input
% Y : n x m is the matrix to be factorized (Y = S M)
% Sinitial = initial value of S
% Minitial = initial value of M

% initializations
[n,m] = size(Y);
r = size(Minitial,1);
p = 2;
k = 0;
epsilon_best = +inf;
sigma = 0.1;
R = zeros(r,m);
gamma = 5;
eta = 0.5;
Kmax = 25;
Lmax = 15;
S = Sinitial;
M = Minitial;



% compute number of frames
f = m/p;

% instrumentation code
% track_obj = [];
% track_epsilon = [];

while k < Kmax
    
    

    
    Ytilde = sqrt(2/sigma)*Y;
    
    l = 0;
    while l < Lmax
        
        % solve subproblem 1
         N = []; L = [];
        for i = 1:f
            [Ni, Li] = proj_deformable_approx(M(:,(i-1)*p+1:i*p)-(1/sigma)*R(:,(i-1)*p+1:i*p));
            N = [ N , Ni ];
            L = [ L , Li ];
        end;        
                  
            
       
        
        % solve subproblem 2    
        C = N + (1/sigma)*R;
        [Q,dummyS,dummyV] = svds([ Ytilde ; C ]',r); 
        
        %Shat = Ytilde*Q; 
        A = C*Q;
        M = A*Q'; 
        %S = Shat/inv(A*sqrt(2/sigma));        
        %fprintf('k = %d     cond(A) = %f     nC = %f     sigma = %f     cond(N) = %f \n',k,cond(A),cond(C),sigma,cond(N));
        
        
        
        % update l
        l = l+1;
        
    end;

% instrumentation code
%    track_obj = [ track_obj , norm(Y-S*N,'fro')^2 ];
%    figure(1);
%    clf;
%    plot(track_obj); 
%    grid on;
        
    epsilon = norm(M-N,'fro')^2;
    
%    track_epsilon = [ track_epsilon , epsilon];
%    figure(2);
%    clf;
%    semilogy(track_epsilon);
%    grid on;
    
%   drawnow;
    
    if epsilon < eta*epsilon_best
        R = R-sigma*(M-N);
        epsilon_best = epsilon;
    else
        sigma = gamma*sigma;
    end;
    
    % update k
    k = k+1;
    
    if rem(k,5) == 0
        fprintf('%2.2f%% done\n',k/Kmax*100);
    end;
    
end;
    

M = N;    
S = Y*pinv(M);



function [Y,L] = proj_deformable_approx(X)

r = size(X,1);
d = r/3;

A = zeros(3,3);
for i = 1:d
    Ai = X((i-1)*3+1:i*3,:);
    A = A + Ai*Ai';
end;

[U,S,V] = svd(A);

Q = U(:,1:2);

G = zeros(2,2);
for i = 1:d
    Ai = X((i-1)*3+1:i*3,:);
    Ti = Q'*Ai;
    gi = [ trace(Ti) ; Ti(2,1)-Ti(1,2) ];
    G = G + gi*gi';
end;

[U1,S1,V1] = svd(G);

G = zeros(2,2);
for i = 1:d
    Ai = X((i-1)*3+1:i*3,:);
    Ti = Q'*Ai;
    gi = [ Ti(1,1)-Ti(2,2) ; Ti(1,2)+Ti(2,1) ];
    G = G + gi*gi';
end;

[U2,S2,V2] = svd(G);

if S1(1,1) > S2(1,1)
    u = U1(:,1);
    R = [ u(1) -u(2) ; u(2) u(1) ];
else
    u = U2(:,1);
    R = [ u(1) u(2) ; u(2) -u(1) ];
end;

Q = Q*R;

Y = [];
L = [];
for i = 1:d
    Ai = X((i-1)*3+1:i*3,:);
    ti = 0.5*trace(Q'*Ai);
    if i == 1 && ti < 0
        ti = -ti;
        Q = -Q;
    end
    L = [ L ; ti];
    Y = [ Y ; ti*Q ];
end;


