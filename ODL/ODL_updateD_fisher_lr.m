function D_tilde = ODL_updateD_fisher_lr(D_tilde, A, B, Y_range, opts)    
    function [cost_fisher, cost] = calc_cost(D)
        cost_fisher = 0.5*opts.lambda2 * FDDL_discriminative(D', Y_range);
        cost = -trace(B*D') + 0.5*trace(A*D'*D) + cost_fisher;
    end
    
    D = D_tilde;
    K = size(D, 2);
    N = size(D, 1);
	Dold = D;
	iter = 0;
    cost_old = 1e6;
    while (iter < opts.max_iter_D)
        iter = iter + 1;
        for i = 1: K
            if(A(i, i) ~= 0)
                a = (A(i, i) * eye(N) + opts.lambda2 * opts.H_tilde) \ (B(:,i) - D*A(:, i) + A(i, i)*D(:,i));
                D(:,i) = a/(max(norm(a,2),1));
            end
        end
        
        %%check cost function convergence
        if 0 == mod(iter, 5)
            [cost_fisher, cost] = calc_cost(D);
            if abs(cost - cost_old) < 1e-6
                break; 
            end
            cost_old = cost;
            if 0
                costs = [costs, cost];
                cost_fishers = [cost_fishers, cost_fisher];
            end
        end
        Dold = D;
    end
    D_tilde = D;
end 