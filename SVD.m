% SVD
fs=33; fs_label = 40; lw = 2;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
format long;


A = [3, 2; 4, 5; 1, 1];
b = [1;4;1];

N = 100;
x = linspace(-2.5,2.5, N);
y = linspace(-1.5,3.5, N);

R = zeros(N);

for i = 1:N
    for j = 1:N
        temp = [x(i); y(j)];
        r = A*temp - b;
        R(i, j) = norm(r);
    end
end

[X, Y] = meshgrid(x, y);

solution = pinv(A)*b

figure
c = contourf(X, Y, transpose(R), 25);
%c.LineStyle = 'none';
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', fs_label);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', fs_label);
%set(c, 'edgecolor', 'none', 'LineWidth', lw);
box on
hold on
colormap turbo
test = colorbar;
colormap turbo
caxis([min(min(R)) max(max(R))]);
test.FontSize = 25;
test.TickLabelInterpreter = 'latex';
% test.Label.String = "$V(x, y)$ [V] ";
% test.Label.FontSize = 30;
% test.Label.Interpreter = 'latex';
p = plot(solution(1),solution(2),'ro', 'markersize', 10, 'MarkerFaceColor', 'r');
t = text(solution(1)+0.1,solution(2), '$x$', 'Interpreter', 'latex', 'fontsize', fs_label, 'Color','red');
