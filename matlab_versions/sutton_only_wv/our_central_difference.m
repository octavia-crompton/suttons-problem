%----------------------------------------------------
%Note: Computes gradients of a data vector using 
%      central differencing. At edges, forward and 
%      backward differences are used
%----------------------------------------------------
function [dsdz]=our_central_differnce(s,dz);
m=length(s);
dsdz(2:m-1)=(s(3:m)-s(1:m-2))/(2*dz);
%---- derivatives at first and last nodes are computed
%     using forward and backward difference in space
dsdz(1)=(s(2)-s(1))/dz;
dsdz(m)=(s(m)-s(m-1))/dz;