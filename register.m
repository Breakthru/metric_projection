function [WT,t] = register(W)
% Remove centroid from a measurement matrix
%
% Author: Marco Paladini (paladini@dcs.qmul.ac.uk)
% Last Modified: 18/08/2009
% License: GPLv2
%
[f2,p]=size(W);
t=zeros(1,f2);
WT=zeros(f2,p);
for i=1:f2
    %for every row, compute the mean and subtract it 
    t(i)=mean(W(i,:));
    WT(i,:) = W(i,:) - t(i);
end

end