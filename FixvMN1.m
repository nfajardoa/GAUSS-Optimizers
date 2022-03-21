% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% FixvMN1: Solves a system of non-linear equations using a modified Newton Method
%
% Usage:   [x,crit] = FixVMN1(x0,@F)
%
%  Input:  @F  := Handle of the vector valued function F(x), whose
%                 zero, F(x1)=0, is to be computed. 
%
%          x0  := k times 1 vector of starting values
% 
%  Output: x1   := k times 1 vector, the approximate solution to F(x1)=0
%          crit :=         1 vector, where
%                                   crit(1)=0 : normal termination 
%                                          =1 : function evaluation failed,
%                                          =2 : no further decrease in function value possible, 
%                                          =3 : maximum number of iterations exceedes
%                                   crit(2)   : termination criterion: maxc(abs(F(x)))
%                                   crit(3)   : the maximum relative change of x between the last two iterations
%                                   crit(4)   : F(x)'F(x)/2
%                                   crit(5)   : number of iterations
%
% Globals: _MNR_Print = true (false) do (not) print messages to the screen
%
%          _MNR_Global = true (false) do (not) use line search
%
%          _MNR_LU = true (false) do (not) use LU factorization to solve for the Newton step
%
%          _MNR_Stop = 1 (0) do (not) use the default stopping criterium of 1.e-7
%
%          _MNR_ParTol = 1 (0) do (not) use the default parameter convergence criterium
%

function [x1, crit] = FixvMN1(x0,F)
  
  _MNR_Global = true;  % Use GetStep in FixvMN            
  _MNR_Print  = true;  % FixvMN print messages in the screen
  _MNR_LU     = true;  % use LU factorization for solution
  _MNR_Stop   = 1 ;    % use default stopping criterium 
  _MNR_ParTol = 1 ;    % use default parameter convergence criterium 
 
  % INITIALIZE
  
  maxit = 500;
  
  if _MNR_Stop == 1  % stopping criterium
    stopc = 1e-5;
  else
    stopc = _MNR_Stop;
  end 

  if _MNR_ParTol==1 % parameter convergence criterium
    pTol = eps^(2/3);
  else
    pTol = _MNR_ParTol;
  end
  
  x1 = x0;
  typf = vfeval(F,x0);
  
  % ITERATIONS
  
  crit = ones(5,1);
  crit(1) = 0;
  critold = 2;

  while not((crit(5) > maxit));      % start iterations
    
    if _MNR_Print == true
      disp(['Step No:', num2str(crit(1))])
      disp(['Convergence criterion ', num2str(crit(2))])
      disp(['Minimization criterion: ', num2str(crit(4)), " decrease in function value? ", num2str(crit(4)<critold)])
    end
      
    df = CDJac(F,x1,size(x1,2));
    
    if isnan(df)
      crit(1)=1;
      disp("Function evaluation failed")
      return
    end
    
    fx = vfeval(F,x1);
      
    if _MNR_Global == true
      dg = fx.'.*df; 
    end

    if _MNR_LU == true
      dx = LSolve(df,-fx); % use the LU factorization
    else
    dx=-inv(df).*fx;
    end
  
    x2 = x1+dx;
    step1 = 1;
    mstep = MinStep(x1,dx,ones(size(x1,1),1),pTol);
    
    while isnan(vfeval(F,(x1+step1*dx))) ~= 0 #ismissing(vfeval(F,(x1+step1*dx)))
      step1=0.75*step1;
    end
    
    if _MNR_Print == true
      disp(["mStep= ", num2str(mstep)]);
      disp(["Step1= ", num2str(step1)]);
    end
      
    if ((mstep <= 1) && (step1 < mstep))
      crit(1) = 2;
      disp("Termination criterion: max(abs(F(x)))")
      return
    end
    
    dx=step1*dx;
    
    if _MNR_Global == true
      step2 = NRStep(x1,dx,dg,pTol,F);
    else
      step2 = 1;
    end

    if _MNR_Print == true
      disp(["Step2= ", num2str(step2)])
    end
       
    x2 = x1+step2*dx;
    crit(2) = max(abs(vfeval(F,x2)));       
    crit(3) = ParTest(x1,step2*dx,ones(size(x1,1),1));
    
    if (crit(2)<stopc);
      crit(1) = 0;
      disp("Normal termination")
      return
    elseif (step1*step2) < mstep 
        crit(1) = 2;
        disp("No further decrease in function value possible")
      return
    end
    
    if _MNR_Global == true     
      critold = crit(4);  
      crit(4) = (vfeval(F,x2).'.*vfeval(F,x2))/2;    
    end
       
    x1 = x2;
    crit(5) = crit(5)+1;
    
    end
    
    if crit(5) >= maxit
       crit(1) = 3;
       disp("Maximum number of iterations exceeded")
       return
    end
    
end