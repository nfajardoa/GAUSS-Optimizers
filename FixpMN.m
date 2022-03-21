% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% FixpMN: Computes the zero of function using a modified Newton method
%
% Usage: x1=FixpMN(x0,bounds,@f)
%
%         where:     x1: the solution to f(x)=0 in the intervall [a,b],
%                    x0: the starting value of the algorithm
%                bounds: a 1x2 vector with bounds(1)=a and bounds(2)=b,
%                    @f: a handle to the procedure that computes and
%                        returns f(x)
%
% Remarks:  1. tol=1.e-10 is the default precision, i.e. the procedure stops
%              if abs(f(x))<tol.
%
%           2. if the procedure is not able to converge in maxit=5000
%              iterations it stops, returns the current value of x 
%              and prints an error message to the default output window.                         
%
%           3. The procedure uses the MATLAB routine "gradient" to
%              compute the slope of f at x.
%                  
%

function [x1,itn] = FixpMN(x0,bounds,f)
 
  tol = 1.e-10;
  maxit = 5000;

  itn = 1;
  x1 = x0;

  while not((abs(f(x1) < tol ) | (itn > maxit)))

    df = gradient(f,x1);
    dx = -f(x1)/df;
    x2 = x1+dx;
    step = 1;
  
    while (x2<bounds(1) | x2>bounds(2) | f(x2)==NaN(1))
  
      step = step/2;
      x2 = x1 + step*dx;
  
    end

    itn = itn+1;
    x1 = x2;
  
    if itn>maxit
    
    error("Maximum number of iterations exceeded.")
    
    end
  
  end

end