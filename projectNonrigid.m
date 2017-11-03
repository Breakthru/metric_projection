function [Rstief,l]=projectNonrigid(M)
% Search for the closest matrix to M that is in the form
% [l1 Ri | l2 Ri | ... | lk Ri]
% where Ri is Stiefel and l1..lk are the deformable shape coefficients
% Optimality is guaranteed by convex reformulation
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
% M: Current estimate of the motion matrix, size 2F by 3K
%
% Output:
%
% Rstief: Camera matrices stacked vertically, size 2*F by 3. Those are exactly Stiefel
% l: Deformation weights matrix, size F by K


Rstief=zeros(size(M,1),3);
StackC=zeros(3*size(M,1),6);
for f= 1:size(M,1)/2
    
A=zeros(6,6);
        K=size(M,2)/3;
        for i=1:K
            m=vec(M(2*f-1:2*f,3*i-2:3*i)'); %CONTROL transpose
            A=A - m*m';
        end
        
StackC(6*f-5:6*f,:) = A;
end
%% Marko Stošić, Marija Dodig and João Xavier code
fsdp = [];
cvx_quiet(true);
% retrieve problem instance from the stack
    A = StackC(1:6,:);

    % solve the instance by a SDP
    cvx_begin
    
        variable Q(6,6) symmetric;
        minimize(trace(Q*A));
        trace(Q(1:3,1:3)) == 1;
        trace(Q(4:6,4:6)) == 1;
        trace(Q(1:3,4:6)) == 0;
        Q == semidefinite(6,6);
        B = Q(1:3,4:6);
        w = [ B(2,3)-B(3,2) ; B(3,1)-B(1,3) ; B(1,2)-B(2,1) ];
        [ eye(3) - Q(1:3,1:3) - Q(4:6,4:6) w ; w' 1 ] == semidefinite(4);
        
    cvx_end
    
    % store optimal value obtained by SDP
    fsdp = [ fsdp , trace(A*Q) ];
    
    % keep solution of the first SDP to initialize Newton's method
        [U,S,V] = svd(Q);
        Uopt = [ U(1:3,1) U(4:6,1) ]*sqrt(2);
fnewton = fsdp(1);

Rstief(1:2,:)=Uopt';

% starts at solution of the first SDP
k = size(M,1)/2;
U = Uopt;
for frame = 2:k
    
    % retrieve problem instance from the stack
    A = StackC((frame-1)*6+1:frame*6,:);
    
    fobj = U(:)'*A*U(:);
    
    while 1
        
        % compute extrinsic and intrinsic gradients at U
        extG = reshape(2*A*U(:),3,2);
        G = extG-U*extG'*U;
    
        % check if U is stationary point
        if norm(G) < 1e-5 fnewton = [ fnewton , U(:)'*A*U(:) ]; break; end;
       
        % compute Newton direction
        c = [ U(2,1)*U(3,2)-U(3,1)*U(2,2) ; -U(1,1)*U(3,2)+U(3,1)*U(1,2);
              U(1,1)*U(2,2)-U(2,1)*U(1,2) ];
        Basis = [ U(:,2) c zeros(3,1) ; -U(:,1) zeros(3,1) c ];
    
        H = zeros(3,3);
        Proj = eye(3)-U*U';
        for i = 1:3
            Bi = reshape(Basis(:,i),3,2);        
            for j = i+1:3
                Bj = reshape(Basis(:,j),3,2);
    
    %            keyboard;
    
                Gamma = 0.5*(Bi*Bj'+Bj*Bi')*U+U*(Bj'*Proj*Bi+Bi'*Proj*Bj)*0.5;
                H(i,j) = 2*Bi(:)'*A*Bj(:)-trace(extG'*Gamma);
            end;
        end;
        H = H + H';
    
        Proj2 = eye(3)-0.5*U*U';
        g = [];
        for i = 1:3
            Bi = reshape(Basis(:,i),3,2);        
            Gamma = Bi*Bi'*U+U*(Bi'*Proj*Bi);
            H(i,i) = 2*Bi(:)'*A*Bi(:)-trace(extG'*Gamma);
    
            g = [ g ; trace(Bi'*Proj2*G) ];
        end;
    
        d = -H\g;
    	if cond(H) > 1e3
            d = -g;
        end;     
        D = reshape(Basis*d,3,2);
    
        % check if D is a descent direction
        if trace(extG'*D) >= 0 D = -G; end;
    
        % pre-computations for the geodesic
        [Q,R] = qr((eye(3)-U*U')*D,0);    
        C = diag(sign(diag(R)));
        Q = Q*C;
        R = C*R;
    
        [V1,D1] = eig([ U'*D -R' ; R zeros(2) ]);
        iV1 = inv(V1);
    
        % perform backtracking search
        f = U(:)'*A*U(:);
        t = 1;
        Ut = real([ U Q ]*V1*expm(t*D1)*iV1*[ eye(2) ; zeros(2) ]);
        while Ut(:)'*A*Ut(:) > f+0.2*t*trace(D'*extG)
            t = t*0.5;
            Ut = real([ U Q ]*V1*expm(t*D1)*iV1*[ eye(2) ; zeros(2) ]);
            if t < 1e-5 fprintf('\n t is too small! gauss newton search failed, using the (slower) SDP'); 
                Ut=restart_sdp(A); break; end;
        end;
            
        U = Ut;       
    
    end;

    Rstief(2*frame-1:2*frame,:) = U';
    
end;        
        
        %% -- end of Marko, JX, Marija code
        l=zeros(k,K);
        for i  = 1:k
        Ri=Rstief(2*i-1:2*i,:);

                    % 24-9-2008, USE the ambiguity to fix l1>0
        l(i,1) =  (1/2) * trace(M(2*i-1:2*i,1:3)'*Ri);
        if l(i,1) < 0 
            Ri=-Ri;
            Rstief(2*i-1:2*i,:)= - Rstief(2*i-1:2*i,:);
            l(i,1)=-l(i,1);
        end
        for j=2:K
            l(i,j)= (1/2) * trace(M(2*i-1:2*i,3*j-2:3*j)'*Ri);
        end
        end
end


function Uopt = restart_sdp(A);
    % solve the instance by a SDP
    cvx_begin
    
        variable Q(6,6) symmetric;
        minimize(trace(Q*A));
        trace(Q(1:3,1:3)) == 1;
        trace(Q(4:6,4:6)) == 1;
        trace(Q(1:3,4:6)) == 0;
        Q == semidefinite(6,6);
        B = Q(1:3,4:6);
        w = [ B(2,3)-B(3,2) ; B(3,1)-B(1,3) ; B(1,2)-B(2,1) ];
        [ eye(3) - Q(1:3,1:3) - Q(4:6,4:6) w ; w' 1 ] == semidefinite(4);
        
    cvx_end
 
 
    % keep solution of the first SDP to initialize Newton's method
        [U,S,V] = svd(Q);
        Uopt = [ U(1:3,1) U(4:6,1) ]*sqrt(2);
end