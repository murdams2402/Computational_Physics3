function x = solve_NonLinear_CG_it(N, dim, X, Diff, Hess, tol, tol_star, it_max)
% This function implements the conjugate gradient algorithm
% Inputs: 
%           tol: tolerance for iterations
%           iterations: maximum number of itertions
% Returns:
%           x: the solution of the the said linear equation


x = rand(N, dim);

% Declaring residual
r = - double(subs(Diff, X, x));
% Declaring direction
d = r;
alpha = 1;

k=1;

x_temp = x;

while norm(r) > tol && k <= it_max
    
    while norm(alpha.*d) > tol_star
        if transpose(d)*double(subs(Hess,X,x_temp))*d >= 0
            alpha = -(transpose(double(subs(Diff,X,x_temp)))*d)/(transpose(d)*(double(subs(Hess,X,x_temp)))*d);
            x_temp = x_temp + alpha.*reshape(d.', [dim, N]);
            % d = double(subs(Diff, X, x_temp));
        else
            x_temp = x_temp - 0.01.*reshape(d.', [dim, N]);
            % d = double(subs(Diff, X, x_temp));
        end
    end
    
    % Update approximative solution
    x = x_temp;
    % Update residual
    r_new = - double(subs(Diff, X, x));
    % Improvement for the step
    beta = r_new'*r_new/(r'*r);
    % Generating new search direction
    d = r_new + beta .* d;
    
    % Updating the residual for second 'while' test
    r = r_new;
    
    k=k+1;
    
end

end