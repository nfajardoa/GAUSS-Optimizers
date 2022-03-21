% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
%
%   vfeval (Vector Function Evaluation)
%
%   Purpose: Evaluate an array of a vector function handle F(x)
%
%   Usage:   A = vfeval(@F,x0);
%
%   Input:  @F := m times 1 cell of anonymous function handles
%                 (Each function have to define the same variables in the same order)
%                 
%           x0 := n times 1 vector of values to be evaluated in @F
%                 (where n is the number of variables defined in @F)
%                 (The order of the variables must be the same of the ones defined in @F)
%         
%   Output:  A := m times 1 vector of values of F(x0)
%

function A = vfeval(F,x0)

  x0 = num2cell(x0);

    for i = 1:numel(F)  
      A{i} = F{i}(x0{:});
    end

  A = cell2mat(A);
  
end