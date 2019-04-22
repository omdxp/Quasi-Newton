function [x, f, n] = Quasi_Newton(fun, x0, nit, eps, method)
x(:, 1) = x0;
rho     = 0.4;
b       = 0.8;
% Premier approximation de H (l'identité)
[~, ~, h] = fun(x(:, 1));
H_old     = eye(length(h));
for i=1:nit
    % Mise à jour de la direction de descente 
    [~, g, ~] = fun(x(:, i));
    dk        = -inv(H_old) * g;
    % Recherche du pas optimal par le critère d'Armijo
    a     = 1;
    x_new = x(:, i) + a * dk;
    f_new = fun(x_new);
    f     = fun(x(:, i));
    while f_new > f + rho * a * g' * dk
        a     = a * b;
        x_new = x(:, i) + a * dk;
        f_new = fun(x_new);
    end
    % Mise à jour de x
    x(:, i+1) = x_new;
    if(norm(dk) < eps)
        break;
    end
    % Mise à jour de H_old
    [~, g_k, ~]   = fun(x(:, i));
    [~, g_k_1, ~] = fun(x(:, i+1));
    yk = g_k_1 - g_k;
    sk = x(:, i+1) - x(:, i);
    if method == 0 % Méthode de BFGS
        H_new = H_old + (yk*yk')/(yk'*sk) - (H_old*(sk)*sk'*H_old)/(sk'*H_old*sk);
    elseif method == 1 % Méthode de DFP
        B_old = inv(H_old);
        B_new = B_old + (sk*sk')/(yk'*sk) - (B_old*(yk)*yk'*B_old)/(yk'*B_old*yk);
        H_new = inv(B_new);
    else
        error('Erreur sur l''utilisation de la fonction');
    end
    H_old = H_new;
end
n = i;
f = fun(x(:, end));