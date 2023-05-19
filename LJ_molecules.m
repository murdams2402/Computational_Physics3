fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;

epsilon = 1;
sigma = 1;

it_max = 1e6;
tol = 1e-3;
tol_star = 1e-3;

dim = 3;
N = 3;

X = sym('x', [dim, N])

[diff_E, Hess_E] = E(X, N, epsilon, sigma);

x = solve_NonLinear_CG_it(N, dim, X, diff_E, Hess_E, tol, tol_star, it_max);

visualize_molecule(x,3,'1','test')

