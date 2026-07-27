function [D_bar, D_tilde, iter] = ODL_updateD_mod2(D_bar, D_tilde, A, B, Y_range, opts)    
    function [cost_fisher, cost_h, cost] = calc_cost(D, K)
        D_bar_tmp = D(:, 1: opts.k_bar);
        D_tilde_tmp = D(:, 1 + opts.k_bar: K);
        cost_fisher = 0.5 * opts.lambda2 * trace(D_tilde_tmp' * opts.H_tilde * D_tilde_tmp);
        cost = -trace(B*D') + 0.5*trace(A*D'*D) + cost_fisher;
        cost_h = 0;
        if ~opts.isJSTSP
            cost_h = 0.5*opts.lambda3 * trace(D_bar_tmp' * opts.H_bar * D_bar_tmp);    
            cost = cost + cost_h;
        end
    end
    
    D = [D_bar D_tilde];
    K = size(D, 2);
    N = size(D_bar, 1);
	Dold = D;
	iter = 0;
    cost_old = 1e6;
    [cost_fishers, cost_hs, sw_trs, sb_trs, costs] = deal([], [], [], [], []);
    while (iter < opts.max_iter_D)
        iter = iter + 1;
        for i = 1: K
            if(A(i, i) ~= 0)
                if i <= opts.k_bar   
                    if 1 == opts.isJSTSP
                        a = 1.0/A(i, i) * (B(:,i) - D*A(:, i)) + D(:,i);
                    else
                        a = (A(i, i) * eye(N) + opts.lambda3 * opts.H_bar) \ (B(:, i) - D * A(:, i) + A(i, i)*D(:,i));
                    end
                else % for discriminative dictionary
                    a = (A(i, i) * eye(N) + opts.lambda2 * opts.H_tilde) \ (B(:,i) - D*A(:, i) + A(i, i)*D(:,i));
                end
                D(:,i) = a/(max(norm(a,2), 1));
            end
        end
        
        %%check cost function convergence
        if 0 == mod(iter, 5)
            [cost_fisher, cost_h, cost] = calc_cost(D, K);
            if abs(cost - cost_old) < 1e-6
                break; 
            end
            cost_old = cost;
            if 0
                costs = [costs, cost];
                cost_fishers = [cost_fishers, cost_fisher];
                cost_hs = [cost_hs, cost_h];
                cost_incoherences = [cost_incoherences, cost_incoherence];
                [sw_tr, sb_tr] = sw_b_tr(D_bar', Y_range);
                [sw_trs, sb_trs] = deal([sw_trs, sw_tr], [sb_trs, sb_tr]);
                figure; 
                plot(cost_fishers(1:end));
                title("cost\_fishers w.r.t iterations");
                
                figure; 
                plot(cost_hs(1:end));
                title("cost\_hs w.r.t iterations");
                
                figure; 
                plot(cost_incoherences(1:end));
                title("cost\_incoherences w.r.t iterations");
                
                figure; 
                plot(sw_trs);
                title("sw\_trs w.r.t iterations");
                
                figure; 
                plot(sb_trs);
                title("sb\_trs w.r.t iterations");
                
                figure; 
                plot(costs(1:end));
                title("costs w.r.t iterations");
                close all;
            end
        end
        Dold = D;
    end
    D_bar = D(:, 1: opts.k_bar);
    D_tilde = D(:, opts.k_bar+1: end);
end 
function [sw_tr, sb_tr] = sw_b_tr(D, Y_range)
    % columbs of D belong to different classes
    C = size(Y_range, 2) - 1;
    sw = 0;
    sb = 0;
    total = 0;
    for i = 1: Y_range(C)
        total = total + D(:, i);
    end
    avg = total / Y_range(C);
    for c = 1: C
        D_c = D(:, Y_range(c) + 1: Y_range(c + 1));
        sum_c = D_c(:, 1);
        for i = 2: size(D_c, 2)
            sum_c = D_c(:, i) + sum_c; 
        end
        avg_c = sum_c / (Y_range(c + 1) - Y_range(c));
        for i = 1: size(D_c, 2)
            sw = sw + (D_c(:, i) - avg_c) * (D_c(:, i) - avg_c)';
        end
        sb = sb + (avg_c - avg)*(avg_c - avg)';
    end
    sw_tr = trace(sw);
    sb_tr = trace(sb);
end