function s=makeS3fxp(l,B)
% Combine basis shapes and deformation weights into a 3D sequence
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% l: Deformation weights matrix, size F by K where F is the # of frames, K the # of basis shapes
% B: Basis Shapes stacked vertically, size 3K by P where P is the number of points
%
% Output:
%
% s: Resulting 3D sequence, 3*F by P

[nF,k]=size(l);
s=zeros(3*nF,size(B,2));

for n=1:nF
    for j=1:k
        s(3*n-2:3*n,:)=s(3*n-2:3*n,:)+l(n,j)*B(3*j-2:3*j,:);
    end
end

end
    