% Example of definitions

% Anonymous functions (Be careful with ther order)
a = @(x, y) x.^2 + y;
b = @(x, y) y.^3 + x;

% Creation of the cell
F = {a, b}

% Definition of evaluation points
x0 = [2, 1]

% Check correct evaluation
vfeval(F,x0)