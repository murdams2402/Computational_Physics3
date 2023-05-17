function [x, it] = solve_NonLinear_CG_it(A, b, tol, tol_star)
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
alpha = 1;

it = zeros(10, 1);
x_th = A\b;
k=1;

x_temp = x;

while norm(r) > tol
    
    while norm(alpha*d) > tol_star
        alpha = - r'*d/(d'*r*d);
        x_temp = x_temp + alpha.*d;
    end
    
    % Update approximative solution
    x = x_temp;
    % Update residual
    r_new = b - A*x;
    % Improvement for the step
    beta = r_new'*r_new/(r'*r);
    % Generating new search direction
    d = r_new + beta .* d;
    
    % Updating the residual for second 'while' test
    r = r_new;
    
    it(k) = norm(x_th - x)/norm(x_th);
    k=k+1;
    
end

end