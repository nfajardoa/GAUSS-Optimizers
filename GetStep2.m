% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% GetStep2:
%
% Purpose:  Find the step size s so that the Newton-Raphson algorithm
%           always moves in the direction of a (local) minimum of (1/2)(f(x)'f(x))
%
% Usage:    sf = GetStep2(x0,dx0,dg,&f);
%
% Input:    x0  := n times 1 vector, the initial point
%           dx0 := n times 1 vector, the Newton direction
%           dg  := 1 times n vector, the gradient of (1/2)(f'f) at x_0
%           @f  := pointer to the function whose zero x solves the system of equations
%
% Output:   sf  := admissible stepsize
%

function sf = GetStep2(x0,dx0,dg,F)

  % Fixed parameters of the algorithm
  smult = 1.0e-4; 
  smin  = 0.1; 
  smax  = 0.5; 
  stol  = 1.e-9;
  
  % Initialize
  s1   = 1.0;
  amat = zeros(2,2);
  bvec = zeros(2,1);
  g0   = (1/2)*(vfeval(F,x0).'.*vfeval(F,x0));
  dgdx = dg*dx0;  % dg(x0)*dx0
  g1   = (1/2)*(vfeval(F,x0+dx0).'.*vfeval(F,x0+dx0));

  % Try the full Newton step s = 1
    if g1 <= g0 + smult*dgdx
      sf = s1;
    else;
      s = -dgdx/(2*(g1-g0-dgdx));
      if s<smin
        s=smin;
      end
      if s>smax 
        s=smax; 
      end
      
      x1 = x0+s*dx0;
      g2 = (1/2)*(vfeval(F,x1).'.*vfeval(F,x1));
    
    end
    
  s2 = s;

  % Reduce s2 further unless g2 < g0 + s2*smult*dgdx */
  while (g2 > (g0 + smult*s2*dgdx))

    amat(1,1) =  1/(s2^2);
    amat(1,2) = -1/(s1^2); 
    amat(2,1) = -s1/(s2^2); 
    amat(2,2) =  s2/(s1^2);
    bvec(1)   =  g2-s2*dgdx-g0; 
    bvec(2)   =  g1-s1*dgdx-g0;
    ab        =  (amat*bvec)/(s2-s1);
    s         =  (-ab(2) + sqrt((ab(2)^2) - 3*ab(1)*dgdx))/(3*ab(1));

    if s < s2*smin 
      s = s2*smin; 
    end
    
    if s > s2*smax
      s = s2*smax; 
    endif

    tol = sqrt((s*dx0)'(s*dx0))/(1+sqrt(x0'x0));
    
    if tol < stol
      retp(-1.0); 
    end
    
    s1 = s2;
    s2 = s;
    g1 = g2;
    x1 = x0+s2*dx0;
    g2 = (1/2)*(f(x1)'f(x1));
    sf = s2;
    
  end
end
