function [x, it] = solve_CG_it(A, b, tol)
% This function implements the conjugate gradient algorithm
% Inputs: 
%           A: input square symmetric positive-definite matrix 
%           b: right-hand side vector of a linear equation A·x = b
%           tol: tolerance for iterations
%           iterations: maximum number of itertions
% Returns:
%           x: the solution of the the said linear equation

x = rand(length(A), 1);

% Declaring residual
r = b - A*x;
% Declaring direction
d = r;

it = zeros(10, 1);
x_th = A\b;
k=1;

while norm(r) > tol && k <= 100
    % Calculating step length
    alpha = r'*r/(d'*A*d);
    % Update approximative solution
    x = x + alpha.*d;
    % Update residual
    r_new = r - alpha.*A*d;
    % Improvement for the step
    beta = r_new'*r_new/(r'*r);
    % Generating new search direction
    d = r_new + beta .* d;
    
    % Updating the residual for 'while' test
    r = r_new;
    
    it(k) = norm(x_th - x)/norm(x_th);
    k=k+1;
    
end

end