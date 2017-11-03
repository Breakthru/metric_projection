function visualise_articulated(W2d,rec3D,axis3D,ind1,ind2)
% Script to visually evaluate the result of a
% 3D reconstruction
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
% Input:
%
% W2d: Input 2D tracking data
% gt3D: 3D ground truth sequence
% rec3D: Recovered deformable 3D sequence

h=figure('Position',[1 1 500 700])


set(gcf, 'color', [1 1 1]);
T = size(rec3D,1)/3;
min_x=min(min(W2d(1:2:end,:)));
max_x=max(max(W2d(1:2:end,:)));
min_y=min(min(W2d(2:2:end,:)));
max_y=max(max(W2d(2:2:end,:)));

mnx=min(min(rec3D(1:3:end,:)));
mny=min(min(rec3D(2:3:end,:)));
mnz=min(min(rec3D(3:3:end,:)));
mxx=max(max(rec3D(1:3:end,:)));
mxy=max(max(rec3D(2:3:end,:)));
mxz=max(max(rec3D(3:3:end,:)));

for t=1:T
    figure(h)
    subplot(2,1,1)
    plot(W2d(2*t-1,:), -W2d(2*t,:), 'g.');
    title('Input 2D tracks')
    axis ([min_x max_x -max_y -min_y])
    %axis equal
    subplot(2,1,2);
    hold off
    plot3(rec3D(3*t-2,ind1), rec3D(3*t-1,ind1), rec3D(3*t,ind1), 'co');
    hold on;
    plot3(rec3D(3*t-2,ind2), rec3D(3*t-1,ind2), rec3D(3*t,ind2), 'mo');
    plot3(axis3D(3*t-2,:),axis3D(3*t-1,:),axis3D(3*t,:),'g-','LineWidth',2);
    %title('3D shape')
    axis ([mnx mxx mny mxy mnz mxz])
    view(0,-90)
    legend('reconstruction (Obj 1)','reconstruction (Obj 2)','reconstruction (Axis)','location','NorthEastOutside');
   drawnow;
   
  
   if 0,
      I = getframe;
      str = sprintf('frame%04d', t);
      imwrite(I.cdata, [str '.jpg'], 'Quality', 100);      
   end
end


end
