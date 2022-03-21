% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% CheckBounds2:
%
% Purpose:  Finds smallest s1 so that x+s1*dx is within bounds 
%
% Usage:    s1 = CheckBounds2(x,dx,bounds)
%
% Input:   x := k times 1 vector, the given value of x0
%         dx := k times 1 vector, the Newton step from x0 to x1
%     bounds := k times 2 matrix, the first (second) colum store the lower (upper) bounds of x
%
% Output:  s := scalar
%

function s1 = CheckBounds2(x,dx,bounds)
  
  l1 = (bounds(:,1)-x)./dx;
  s1 = min(l1(l1(:,:)>0));
  l2 = (bounds(:,2)-x)./dx;
  s2 = min(l2(l2(:,:)>0));
  s1 = min([0.98*s1,0.98*s2,1]);
 
end