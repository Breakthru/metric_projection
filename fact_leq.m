% Created by Alessio Del Bue and Lourdes Agapito
% Last Modified: 18/09/2003
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
