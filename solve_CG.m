function x = solve_CG(A, b, tol, iterations)
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

k = 1;
while norm(r) > tol && k <= iterations
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
    
    k = k + 1;
end

end