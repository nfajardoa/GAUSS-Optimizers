% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
%@ ----------------------------------- ParTest -------------------------------------------------- 
%
%   Purpose: Compute relative change in x
%
%   Usage:   crit=ParTest(x,dx,typx);
%
%   Input:   x := n times 1 vector
%           dx := n times 1 vector, change in x
%         typx := n times 1 vector, typical elements of x
%
%   Output: crit := scalar
%
%----------------------------------------------------------------------------------------------- @

function maxCrit = ParTest(x,dx,typx)
  
  n=size(x,1);
  crit = zeros(n,1);
  i = 1;
  
  while not(i>n)
    
    crit(i) = abs(dx(i))/max([abs(x(i)),abs(typx(1))]);
    i = i+1;
 
  end
  
 maxCrit = max(crit);

end