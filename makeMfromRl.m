function M = makeMfromRl(l,Rstief)
% This function combines 2 by 3 camera rotations with deformation
% weights to produce the 2F by 3K  Motion matrix
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% l, Rstief: weights and rotations
% Output:
%
% M: Motion matrix

[f,K]=size(l);f2=2*f;
M=zeros(f2,3*K);
    for f=1:2:f2
	for j=1:K 
        M(f:f+1,3*j-2:3*j)=l(ceil(f/2),j)*Rstief(f:f+1,:); 
    end
    end
    
end