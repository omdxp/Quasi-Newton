function [x, f, n] = Quasi_Newton(fun, x0, nit, eps, method)
x(:, 1) = x0;
rho     = 0.4;
b       = 0.8;
% First approximation for H matrix
[~, ~, h] = fun(x(:, 1));
H_old     = eye(length(h));
for i=1:nit
    % Update descent direction
    [~, g, ~] = fun(x(:, i));
    dk        = -inv(H_old) * g;
    % Searching for optimal step by Armijo's method
    a     = 1;
    x_new = x(:, i) + a * dk;
    f_new = fun(x_new);
    f     = fun(x(:, i));
    while f_new > f + rho * a * g' * dk
        a     = a * b;
        x_new = x(:, i) + a * dk;
        f_new = fun(x_new);
    end
    % Update x
    x(:, i+1) = x_new;
    if(norm(dk) < eps)
        break;
    end
    % Update H_old
    [~, g_k, ~]   = fun(x(:, i));
    [~, g_k_1, ~] = fun(x(:, i+1));
    yk = g_k_1 - g_k;
    sk = x(:, i+1) - x(:, i);
    if method == 0 % BFGS method
        H_new = H_old + (yk*yk')/(yk'*sk) - (H_old*(sk)*sk'*H_old)/(sk'*H_old*sk);
    elseif method == 1 % DFP method
        B_old = inv(H_old);
        B_new = B_old + (sk*sk')/(yk'*sk) - (B_old*(yk)*yk'*B_old)/(yk'*B_old*yk);
        H_new = inv(B_new);
    else
        error('Error using function');
    end
    H_old = H_new;
end
n = i;
f = fun(x(:, end));
