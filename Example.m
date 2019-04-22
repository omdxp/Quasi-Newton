clc, clear, close all;
%% 1) Méthode de Newton
% La fonction rosenbrock:
% [x,y] = meshgrid(-1.5:.1:1.5, -1.5:.1:1.5);
% f = (1 - x).^2 + 10 * (y - x.^2).^2;
% v = [0.01,0.5,2,5,10,15,25,35,50,100,200];
% [C, h] = contour(x, y, f, v);
% clabel(C, h);
% grid; hold on;
% % Cholosky method for Ax = b
% % A+ -> p=0
% % b)
% nit = 100;
% x0 = [0; 15];
% eps = 10e-7;
% tic
% [x, f, n] = Newton_gradient(@fun_obj, x0, nit, eps);
% toc
% plot(x(1, :), x(2, :), 'm-d',...
%     'LineWidth', 1.5,...
%     'MarkerEdgeColor', 'b',...
%     'MarkerFaceColor', 'b');
%% 2) Méthode de Quasi-Newton
figure;
clear x;
% La fonction rosenbrock:
[x,y] = meshgrid(-1.5:.1:1.5, -1.5:.1:1.5);
f = (1 - x).^2 + 10 * (y - x.^2).^2;
v = [0.01,0.5,2,5,10,15,25,35,50,100,200];
[C, h] = contour(x, y, f, v);
clabel(C, h);
grid; hold on;
nit = 100;
x0 = [0; 15];
eps = 10e-7;
% 0 pour BFGS, 1 pour DFP
tic
[x, f, n] = Quasi_Newton(@fun_obj, x0, nit, eps, 1);
toc
plot(x(1, :), x(2, :), 'm-d',...
    'LineWidth', 1.5,...
    'MarkerEdgeColor', 'b',...
    'MarkerFaceColor', 'b');