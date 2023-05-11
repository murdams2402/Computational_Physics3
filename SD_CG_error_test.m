fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;


data = load('Matrices.mat');

% data.A1
N = 10;
tol = eps;

x_th = data.A1\data.b1;
err_SD = zeros(N, 1);
err_CG = zeros(N, 1);
N_itr = 10.^(1:N);

for i = 1:N
    x_SD = solve_SD(data.A1, data.b1, tol, 10^i);
    x_CG = solve_CG(data.A1, data.b1, tol, 10^i);
    err_SD(i) = norm(x_th - x_SD)/norm(x_th);
    err_CG(i) = norm(x_th - x_CG)/norm(x_th);
end


figure
plot(N_itr, err_SD, '-+b', 'Linewidth', lw)
xlabel('$x$ [nm]', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$ [nm]', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on
hold on
plot(N_itr, err_CG, '-+m', 'Linewidth', lw)


