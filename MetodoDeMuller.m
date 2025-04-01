function raiz = MetodoDeMuller()
    f = @(x) x.^3 - 13*x - 12; % Definir la función
    x0 = 4.5; 
    x1 = 5.5; 
    x2 = 5; 
    tol = 1e-6; 
    max_iter = 100;

    raiz = muller(f, x0, x1, x2, tol, max_iter);
end

function raiz = muller(f, x0, x1, x2, tol, max_iter)
    for iter = 1:max_iter
        h1 = x1 - x0;
        h2 = x2 - x1;
        delta1 = (f(x1) - f(x0)) / h1;
        delta2 = (f(x2) - f(x1)) / h2;
        d = (delta2 - delta1) / (h2 + h1);

        b = delta2 + h2 * d;
        D = sqrt(b^2 - 4*f(x2)*d);

        if abs(b - D) < abs(b + D)
            E = b + D;
        else
            E = b - D;
        end

        h = -2 * f(x2) / E;
        x3 = x2 + h;

        fprintf('Iteración %d: x = %.6f\n', iter, x3);

        if abs(h) < tol
            raiz = x3;
            return;
        end

        x0 = x1;
        x1 = x2;
        x2 = x3;
    end

    error('El método de Müller no convergió en %d iteraciones', max_iter);
end
