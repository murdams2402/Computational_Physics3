function x = solve_NonLinear_CG_it(x_init, N, dim, X, Diff_E, Hess_E, tol, tol_star)
% This function implements the conjugate gradient algorithm
% Inputs: 
%           x_init: initial guess for the minimum
%           N: number of atoms
%           dim: dimension of the space
%           X: N x dim sym object matrix, which contains the vector 
%              coordinates of each atom in it's rows
%           Diff_E: 1 x N·dim sym object column vector, gradient of the
%           energy function
%           Hess_E: N·dim x N·dim sym object, Hessian matrix of the energy
%           function
%           tol: tolerance for iterations
%           tol_star: second tolerance for inner iterations
% Returns:
%           x: N x dim matrix that minimizes E


x = x_init; %rand(N, dim);

% Declaring residual
r = - double(subs(Diff_E, X, x));
r_new = r;
% Declaring direction
d = r;
while norm(r) > tol
    alpha = -(double(subs(Diff_E,X,x))'*d)/(d'*double(subs(Hess_E,X,x))*d);
    while norm(alpha*d) > tol_star
        if d'*double(subs(Hess_E,X,x))*d >= 0
            alpha = -( double(subs(Diff_E,X,x)'*d)/(d'*(double(subs(Hess_E,X,x)))*d) );
            x = x + alpha*(reshape(d.', [dim, N])).';
            d = double(subs(Diff_E, X, x));
        else
            x = x - 0.01*(reshape(d.', [dim, N])).';
            d = double(subs(Diff_E, X, x));
        end
    end
    % Update residual
    r_new = - double(subs(Diff_E, X, x));
    % Improvement for the step
    beta = r_new'*r_new/(r'*r);
    % Generating new search direction
    d = r_new + beta .* d;
    
    % Updating the residual
    r = r_new;
end

end