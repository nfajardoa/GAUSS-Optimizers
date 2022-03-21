% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% FixvMN: Solves a system of non-linear equations using a modified Newton Method
%
% Usage:   [fx1, x1] = FixvMN(x0,bounds,@F)
%
%  Input:  @F  := Handle the vector valued function F(x), whose
%                 zero, F(x1) = 0, is to be computed. 
%
%          x0  := k times 1 vector of starting values
% 
%      bounds  := k times 2 matrix, whose first column holds the lower and whose
%                                   second column holds the upper bounds on x
%
%  Output: fx1 := k times 1 vector, the value of F at the solution x1
%
%          x1  := k times 1 vector, the approximate solution to F(x1) = 0
%

function [fx1, x1] = FixvMN(x0,bounds,F)

  _MNR_Global = false;  % Do not use GetStep in FixvMN             
  _MNR_Print  = true;   % FixvMN print messages to the screen   
  _MNR_QR     = false;  % do not use QR factorization for solution 

  % Initialize */
  maxit = 5000;     % stop after 5000 Iterations
  stopc = 1e-8;     % stopping criterium
  x1    = x0;

  % Start Iterations
  itn   = 1;
  crit  = 1;
  crit1 = 1;
  crit2 = 1;
  
  while not(((crit<stopc) | (itn > maxit)))      % start iterations
    
    if _MNR_Print == true
      
      disp(["Step No: ", num2str(itn), " Convergence criterion: ", num2str(crit)]);
      
      if _MNR_Global == true
        disp(["Minimization criterion: ", num2str(crit1), " crit1<crit2: ", num2str(crit1<crit2)]);
      end
    
    end
    
    df = CDJac(F,x1,size(x1,2));
    
    if isnan(df)
      disp("Not able to compute gradient");
      return
    end
    
    if _MNR_Global == true
      dg = vfeval(F,x1).'.*df; 
    end

    if _MNR_QR == true
      dx = LSolve(df,-f(x1)); % use the QR-factorization
      if isnan(dx)
        disp("Not able to solve for dx")
        return
      end
    else
      dx = -inv(df).*vfeval(F,x1);  
    end    
    
    x2    = x1 + dx;        
    step1 = CheckBounds2(x1,dx,bounds); % scale down by step if array bounds are exceeded
      
    if _MNR_Print == true
      disp(["Step1= ", num2str(step1)]);
    end
      
    if step1<1.0e-16;
      disp("Not able to find x within admissible bounds")
      disp("Results may be missleading")
      return
    end
    
    dx = step1*dx;
    
    if _MNR_Global == true
      step2 = GetStep2(x1,dx,dg,F);
    else 
      step2 = step1; 
    end
      
    if _MNR_Print == true
      disp(["Step2= ", num2str(step2)]);
    end
    
    if step2<0.0
      disp("Stepsize cannot be further reduced."); 
      disp("Algorythm may stuck at a local minimum");
      disp("Results may be misleading");
      return;
    end
    
    x2   = x1+step2*dx;
    crit = max(abs(vfeval(F,x2)));
    
    if _MNR_Global == true
      crit2 = crit1;
      crit1 = (vfeval(F,x2).'.*vfeval(F,x2))/2;    
    end;
    
    x1  = x2;
    fx1 = vfeval(F,x1);
    itn = itn+1;
   
   end
   
   if itn >= maxit
    disp("Maximum number of iterations exceeded! Returned solution may be missleading");
    return
   end

 end