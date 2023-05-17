fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;

% N has to be an odd number
N = 45;
L = 45;

[H, b] = generate_Hb(N);

figure
s = surf(H'*H);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap gray

x = H\b;

V = zeros(N, N);

for i = 1:N
    for j= 1:N
        V(i, j) = x(index(i, j, N));
    end
end

h=L/(N-1);
[X,Y] = meshgrid(0:h:L);

figure
s = surf(X, Y, transpose(V));
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V)) max(max(V))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';


%% Testing with SD and CG methods
[x_SD, err_SD] = solve_SD_it(H'*H, H'*b, eps);
%%
[x_CG, err_CG] = solve_CG_it(H'*H, H'*b, eps);
%%

V_SD = zeros(N, N);
V_CG = V_SD;

for i = 1:N
    for j= 1:N
        V_SD(i, j) = x_SD(index(i, j, N));
        V_CG(i, j) = x_CG(index(i, j, N));
    end
end

figure
s = surf(X, Y, transpose(V_SD));
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V_SD)) max(max(V_SD))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';


figure
s = surf(X, Y, transpose(V_CG));
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V_CG)) max(max(V_CG))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';


n_SD = length(err_SD);
n_CG = length(err_CG);

N_SD = (1:n_SD);
N_CG = (1:n_CG);

figure('Name', 'SD')
semilogy(N_SD, err_SD, '-+b', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on

figure('Name', 'CG')
semilogy(N_CG, err_CG, '-+m', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on

%% Plotting the electric field (using exact solution for the potential V)
x = (0:h:L);
[Ex, Ey] = gradient(V);
% [Ex, Ey] = gradient(V_CG) or gradient(V_SD)
Ex = -Ex; Ey=-Ey;

scale = 3.5;

figure
% contour, contourfm pcolor or surf
c = pcolor(X, Y, transpose(V));
%c.LineStyle = 'none';
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(c, 'edgecolor', 'none', 'LineWidth', lw);
box on
hold on
q=quiver(x, x, Ex, Ey, scale, '-k');
q.LineWidth = lw;
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V)) max(max(V))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';

%% Jacobi preconditioning
A = H'*H;
n = length(A);
M = zeros(n, n);
for i = 1:n
    M(i, i) = A(i, i);
end
M_inv = M^-1;

[x_SD, err_SD] = solve_SD_it(M_inv*H'*H, M_inv*H'*b, eps);
[x_CG, err_CG] = solve_CG_it(M_inv*H'*H, M_inv*H'*b, eps);


V_SD = zeros(N, N);
V_CG = V_SD;

for i = 1:N
    for j= 1:N
        V_SD(i, j) = x_SD(index(i, j, N));
        V_CG(i, j) = x_CG(index(i, j, N));
    end
end

figure
s = surf(X, Y, transpose(V_SD));
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V_SD)) max(max(V_SD))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';


figure
s = surf(X, Y, transpose(V_CG));
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', fs_label);
set(s, 'edgecolor', 'none', 'LineWidth', lw);
box on
colormap jet
test = colorbar;
colormap jet
caxis([min(min(V_CG)) max(max(V_CG))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
test.Label.String = "$V(x, y)$ [V] ";
test.Label.FontSize = 30;
test.Label.Interpreter = 'latex';


n_SD = length(err_SD);
n_CG = length(err_CG);

N_SD = (1:n_SD);
N_CG = (1:n_CG);

figure('Name', 'SD')
semilogy(N_SD, err_SD, '-+b', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on

figure('Name', 'CG')
semilogy(N_CG, err_CG, '-+m', 'Linewidth', lw)
xlabel('$N_{\rm iterations}$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$\epsilon$', 'Interpreter', 'latex', 'fontsize', fs_label);
box on
grid on



