function [Q,scale,a]=align_shapes(s1,s2,plot)
% Compare two 3D shapes
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% s1: test shape
% s2: reference shape
% plot: option to plot comparisons, 0 - none
%				    1 - before scaling
%				    2 - after scaling
%
% Output:
%
% Q: Rotation that align s1 to s2
% scale: Scale factor, point distances in s1/scale should be equal to those in s2
% a: Shape result after scaling and alignement

if size(s1,1)~=3 || size(s2,1)~=3
    error('align_shapes, shapes size is not 3xP.')
end
if nargin<3 plot=0; end

p1=size(s1,2);
p2=size(s2,2);
if p1~=p2
    error('align_shapes, the number of points is not the same.')
end

s1=register(s1);
s2=register(s2);
if plot~=0
    h=figure()
end

if nargin<3
    plot=0;
end
% plottem
if plot==2
    h=figure()
    plotshape(s1,'b*');
    hold on
    plotshape(s2,'r*');
end

%points between -1 and 1
%s1=s1/max(max(abs(s1)));

%calculate the scale based on point distances
d1=[]; d2=[];
%for p = [1:p1-1]  %mean distance between all points
for p = 1
    d1(p)=norm(s1(:,p)-s1(:,p+1));
    d2(p)=norm(s2(:,p)-s2(:,p+1));
end

scale=d1/d2;
 s1=s1/scale;



if (plot==2)
    plotshape(s1,'b+');
        hold on
    plotshape(s2,'r+');
end

% check if they superimpose
[a,Q,val]=procrust(s2',s1');
Q=Q';
a=a';

if (plot==1)
    hold off
    plotshape(a,'g*');
        hold on
    plotshape(s2,'ko');
    legend('Recovered', 'Ground Truth','Location','Best')
    cameratoolbar
end


%end of function
end

function plotshape(shape,linespec)
    plot3(shape(1,:),shape(2,:),shape(3,:),linespec);
end