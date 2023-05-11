function [x, it] = solve_SD_it(A, b, tol)
% This function implements the steepest descent algorithm
% Inputs: 
%           A: input square matrix 
%           b: right-hand side vector of a linear equation A·x = b
%           tol: tolerance for iterations
%           iterations: maximum number of itertions
% Returns:
%           x: the solution of the the said linear equation

x = rand(length(A), 1);

% Declaring residual
r = b - A*x;

it = [];
x_th = A\b;

k=1;
while norm(r) > tol 
    % Calculating step length
    alpha = transpose(r)*r/(transpose(r)*A*r);
    % Update approximative solution
    x = x + alpha.*r;
    % Update residual
    r = b - A*x;
    
    it(k) = norm(x_th - x)/norm(x_th);
    k=k+1;
end

end