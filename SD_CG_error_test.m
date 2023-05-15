fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;


data = load('Matrices.mat');

%% A1 matrix
tol = eps;

A = data.A1;
b = data.b1;

[~, err_SD] = solve_SD_it(A, b, tol);
[~, err_CG] = solve_CG_it(A, b, tol);


N_itr = (1:length(err_SD));

figure
semilogy(N_itr, err_SD, '-+b', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on

N_itr = (1:length(err_CG));
figure
semilogy(N_itr, err_CG, '-+m', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on


%% Matrix A2

tol = eps;

A = data.A2;
b = data.b2;

[~, err_SD] = solve_SD_it(A, b, tol);
[~, err_CG] = solve_CG_it(A, b, tol);


N_itr = (1:length(err_SD));

figure
semilogy(N_itr, err_SD, '-+b', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on

N_itr = (1:length(err_CG));
figure
semilogy(N_itr, err_CG, '-+m', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on
%%
 kappa = cond(A)
% SD
n_eps = log(eps)/log((kappa-1)/(kappa+1))


