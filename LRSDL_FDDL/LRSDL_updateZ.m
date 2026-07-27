function [Z_bar, Z_tilde] = LRSDL_updateZ(Y, D_bar, D_tilde, Z_bar, Z_tilde, opts)
    function cost = calc_f(Z)
        [Z_bar, Z_tilde] = extractFromZ(Z);
        % fprintf("%f %f %f %f\n", size(D_tilde, 1), size(D_tilde, 2), size(Z_tilde, 1), size(Z_tilde, 2));
        Ybar = Y - D_tilde * Z_tilde;
        cost = 0.5*normF2(Ybar - D_bar * Z_bar);
    end

    function cost = calc_F(Z)
        cost = calc_f(Z) + opts.lambda1*norm1(Z);
    end

    function [Z_bar, Z_tilde] = extractFromZ(Z)
        Z_bar = Z(1: opts.k_bar, :);
        Z_tilde = Z(opts.k_bar + 1: end, :);
    end    

    %% Gradient for FISTA 
    function g1 = grad(Z)
        [Z_bar, Z_tilde] = extractFromZ(Z);
        g1 = -DtY + DtD*Z;
    end
    D = [D_bar D_tilde];
    DtY = D'*Y;
    DtD = D'*D;
    %% check gradient
    if opts.check_grad && ~check_grad(@calc_f, @grad, [Z_bar; Z_tilde])
        fprintf('Check gradient or cost again!\n')
        pause
    end       
    %% ========= Main FISTA ==============================
    L = max(eig(DtD));
    Z = [Z_bar; Z_tilde];

    Z = fista(@grad, Z, L, opts.lambda1, opts, @calc_F);
    [Z_bar, Z_tilde] = extractFromZ(Z);
end 

