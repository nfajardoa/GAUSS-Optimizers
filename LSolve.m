% Written for MATLAB by Nicolas Fajardo & Eduardo Fajardo --- Dec 2018
% Taken from Dynamic General Equilibrium Modelling by Heer, B. & Maussner, A.
%
% LSolve:
%
%  Purpose: Solve a system of linear equations Ax=b using the QR-decomposition
%
%  Usage:  x=LSolve(A,b)
%
%  Input:  A: square matrix of dimension n
%          b: vector of dimension n
%
%  Output: x: vector of dimension n
%   

function x = LSolve(A,b)
  
  %uptriang = dsp.UpperTriangularSolver
  
  [q,r] = qr(A);
  %x = uptriang(q.'.*b,(q.'.*q).*r);
  opts.UT = true;
  x = linsolve(q.'.*b,(q.'.*q).*r);
  
end