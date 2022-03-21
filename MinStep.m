% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
%--------------------------------------- MinStep ----------------------------------------------------
%
%  Purpose: Compute the minimal step size so that x1=x0+mstep*dx=0
%           (in terms of the parameter tolerance criterium pTol)
%
%  Usage: mstep=MinStep(x,dx,typx,pTol);
%
%  Input:  x := n times 1 vector
%         dx := n times 1 vector, change in x
%        typx:= n times 1 vector, typical elements of x          
%        pTol:= scalar
%
%  Output: mstep:= scalar
%---------------------------------------------------------------------------------------------------   @

function Step = MinStep(x,dx,typx,pTol)
  
  n = size(x,1);
  i = 1;
  converge = 0;
  
  while not(i>n)
    
    temp = abs(dx(i))/max([abs(x(i)),abs(typx(1))]);
    if (temp > converge)
      converge = temp;
    end
    
    i = i+1;
    
  end
  
  Step = pTol/converge;
  
end