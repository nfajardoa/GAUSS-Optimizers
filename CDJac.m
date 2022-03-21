% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% CDJac
%
% Usage: df = CDJac(@F,x0,n)
%
%
% Purpose: computes a central difference approximation of the Jacobian
%          matrix of a system of n non-linear functions y_i=f^i(x), where
%          x is a column vector of dimension m.
%
% Input:  @F := handle to the routine that returns the m-vector f(x)
%         x0 := m vector x0, the point at which the derivatives are to be evaluated
%         n  := the size of the vector F(x)
%
% Output: dx := n by m matrix of partial derivatives
%
% Algorithm: based on (A.2.8) in Heer and Maussner, see also Dennis and Schnabel (1983), 
%            Algorithm A5.6.4.
%

function df = CDJac(F,x0,n)
  
  m = size(x0,2);
  df = zeros(m,n);
  Macheps = eps^(1/3);
  x1 = x0;
  x2 = x0;
  i = 1;

  while not(i>n)
    
    if x0(i)<0
      h = -Macheps*max([abs(x0(i)),1]);
    else
      h =  Macheps*max([abs(x0(i)),1]);
    end
    
  temp  = x0(i);
  x1(i) = temp+h;
  x2(i) = temp-h;
  h  = x1(i)-temp; % Trick to increase precision slightly, see Dennis and Schnabel (1983), p. 99 @
  f1 = vfeval(F,x1);
  f2 = vfeval(F,x2);
  j  = 1;
  
    while not(j>m)
      df(j,i) = (f1(j)-f2(j))/(2*h);
      j = j+1;
    end
  
  x1(i)=x0(i);
  x2(i)=x0(i);
  i=i+1;
  
  end
end