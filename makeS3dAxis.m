function [ourS3D, axis3D,Wrepr] = makeS3dAxis(M_,S_,t1,t2,ind1,ind2)
% Articulated Hinge joint recovery by Factorization
% Reconstruction function based on Tresadern-Reid method for hinge
% joint recovery, see reference for further details
%
% Author: Marco Paladini, João Fayad, Alessio Del Bue, Lourdes Agapito
% Last Modified: 6/10/2009
%
% Ref: "Articulated Structure from Motion by Factorization"
% Tresadern, P. and Reid, I.
% CVPR '05: Proceedings of the 2005 IEEE Computer Society Conference on Computer Vision and Pattern Recognition
% 
% Input:
% M_: Resulting Motion matrix from articulated factorization
% S_: Combined Shape matrix of the two objects, rank 5
% t1: Centroid coordinates for first object
% t2: Centroid coordinates for second object
% ind1: Array of the indexes for the points in W belonging to the first object
% ind2: Array of the indexes for the points in W belonging to the second object
%
% Output:
% ourS3D: Reconstructed 3D shape of both objects.
% axis3D: 3D coordinates of points lying on the rotation axis.
% Wrepr: reconstructed 2d tracks of feature points, each frame stacked vertically, each
% 	column represents a point



RRtt = [M_,t2'-t1'];
[A,B,C] = svd(RRtt);
dd_1 = C(:,end);
dd_1 = dd_1 / -dd_1(end);
u1plusu2=dd_1(1);
v1 = dd_1(2);
w1 = dd_1(3);
v2 = dd_1(4);
w2 = dd_1(5);
l1=t1'+M_(:,1:3)*[0;v1;w1];
l2=t1'+M_(:,1:3)*[0.5;v1;w1];
R1=M_(:,1:3);
R2=M_(:,[1 4 5]);
S1 = S_(1:3,ind1);
S2 = S_([1 4 5],ind2);
F=size(M_,1);
ourS3D=[];
axis3D=[];
axisptsX = [-250:10:250];
Saxis = [axisptsX;zeros(1,size(axisptsX,2))+v1;zeros(1,size(axisptsX,2))+w1];
for f=1:2:F
r31=cross(M_(f,1:3),M_(f+1,1:3));
r32=cross(M_(f,[1 4 5]),M_(f+1,[1 4 5]));
% rot angle matrix = inv [c2 c3] * [c2' c3']
%anglem = inv(M_(f:f+1,2:3))*M_(f:f+1,4:5)
S1f=[R1(f:f+1,:);r31]*S1 + repmat([t1(f);t1(f+1);1],1,size(ind1,2));
%     plot3(S(1,:),S(2,:),S(3,:),'r*'); hold on; % first shape
S2f=[R2(f:f+1,:);r32]*S2 + repmat([R1(f:f+1,:);r31]*[u1plusu2;v1;w1],1,size(ind2,2)) + repmat([R2(f:f+1,:);r32]*[0;v2;w2],1,size(ind2,2)) + repmat([t1(f);t1(f+1);1],1,size(ind2,2));
%      plot3(S(1,:),S(2,:),S(3,:),'b*'); % second shape
S=[R1(f:f+1,:);r31]*Saxis + repmat([t1(f);t1(f+1);1],1,size(Saxis,2));
%      plot3(S(1,:),S(2,:),S(3,:),'k-');
  ourS3D = [ourS3D; [S1f S2f]];
  axis3D = [axis3D; S];
end


% compute repr err
S1 = S_(1:3,ind1);
S21 = S_(1,ind2) + u1plusu2;
Sv = repmat(v1,size(S21));
Sw = repmat(w1,size(S21));
S22 = [S_(4,ind2) + v2;S_(5,ind2) + w2];
Z = zeros(2,size(ind1,2));
O = ones(1,size(S_,2));
C1 = [S1;Z];
C2 = [S21;Sv;Sw;S22];

Srepr = [C1,C2;O];
Wrepr= [M_ t1']*Srepr;

end