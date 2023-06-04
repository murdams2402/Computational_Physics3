function [x, it] = solve_CG_JacobiPrecond_it(A, b, tol)
% This function implements the steepest descent algorithm
% Inputs: 
%           A: input square symmetric positive-definite matrix 
%           b: right-hand side vector of a linear equation A·x = b
%           tol: tolerance for iterations
% Returns:
%           x: the solution of the the said linear equation
%           it: error of the output solution x (with respect to the
%           analitycal solution given by Matlab A\b) at each iteration 

x = rand(length(A), 1);

% Jacobi preconditioner
M = diag(diag(A));
M_inv = M^-1;

% Declaring residuals
r = b - A*x;
d = M_inv * r;

% Prealocating space  
it = zeros(100, 1);
x_th = A\b;

k=1;
while norm(r) > tol
    % Calculating step length
    alpha = transpose(r)*M_inv*r/(transpose(d)*A*d);
    % Update approximative solution
    x = x + alpha.*d;
    % Calculating new residuals
    r_new = r - alpha*A*d;
    beta = transpose(r_new)*M_inv*r_new/(transpose(r)*M_inv*r);
    d = M_inv*r_new + beta*d;
    
    % Updating the residual r
    r = r_new;
    
    it(k) = norm(x_th - x)/norm(x_th);
    k=k+1;
end

end