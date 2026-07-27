function [D_bar, D_tilde, iter] = updateD_mod2_nonFinsher(D_bar, D_tilde, A, B, opts)    
    D = [D_bar D_tilde];
    K = size(D, 2);
    N = size(D_bar, 1);
	Dold = D;
	iter = 0;
    cost_old = 1e6;
    while (iter < opts.max_iter_D)
        iter = iter + 1;
        for i = 1: K
            if(A(i, i) ~= 0)
                a = 1.0/A(i, i) * (B(:,i) - D*A(:, i)) + D(:,i);
                D(:,i) = a/(max(norm(a,2),1));
            end
        end
        
        %%check cost function convergence
        if 0 == mod(iter, 5)
            cost = -trace(B*D') + 0.5*trace(A*D'*D);
            if abs(cost - cost_old) < 1e-6
                break; 
            end
            cost_old = cost;
        end
        Dold = D;
    end
    D_bar = D(:, 1: opts.k_bar);
    D_tilde = D(:, opts.k_bar+1: end);
end 