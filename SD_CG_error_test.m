fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;


data = load('Matrices.mat');

% data.A1
N = 10;
tol = eps;

A = data.A1;
% A = data.A2;
b = data.b1;
% b = data.b2;

[~, err_SD] = solve_SD_it(A, b, tol);
[~, err_CG] = solve_CG_it(A, b, tol);


figure
plot(N_itr, err_SD, '-+b', 'Linewidth', lw)
xlabel('$x$ [nm]', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$ [nm]', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on
plot(N_itr, err_CG, '-+m', 'Linewidth', lw)


