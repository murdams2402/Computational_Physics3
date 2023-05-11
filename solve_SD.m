function x = solve_SD(A, b, tol, iterations)
% This function implements the steepest descent algorithm
% Inputs: 
%           A: input square matrix 
%           b: right-hand side vector of a linear equation A·x = b
%           tol: tolerance for iterations
% Returns:
%           x: the solution of the the said linear equation

x = rand(length(A), 1);

% Declaring residual
r = b - A*x;
k=1;
while norm(r) > tol && k <= iterations
    % Calculating step length
    alpha = transpose(r)*r/(transpose(r)*A*r);
    % Update approximative solution
    x = x + alpha.*r;
    % Update residual
    r = b - A*x;
    
    k=k+1;
end

end