function [D_bar, D_tilde, Z_bar, Z_tilde, cost_diff] = Solve_permutation_JSTSP(D_bar, D_tilde, Z_bar, Z_tilde, Y_range, opts)
    D = [D_bar, D_tilde];
    Z = [Z_bar; Z_tilde];
    [K_bar, K_tilde] = deal(opts.k_bar, opts.k_tilde); 
    K = K_bar + K_tilde;

    [cost_tilde, ~, ~, ~] = fisher_cost(D_tilde', Y_range);
    cost_old = 0.5 * opts.lambda2 * cost_tilde;

    cost_assigned = zeros(K, K);
    for k = 1: K
        cost_assigned(k, :) = cost_k(D, Y_range, K_bar, K_tilde, k, opts);
    end

    assign = matchpairs(cost_assigned, 1000000);
    % assign(i, 1) is the new position of the original index assign(i, 2)
    
    % permutate the atoms and maps
    D = D(:, assign(:, 1));
    Z = Z';
    Z= Z(:, assign(:, 1))';
    %{
    % another way to permutate the atoms and maps
    permutation = zeros(K, K);
    for k = 1: K
        permutation(assign(k), k) = 1; 
    end
    D = D * permutation;
    Z = permutation' * Z;
    %}
    [D_bar, D_tilde] = deal(D(:, 1: K_bar), D(:, 1 + K_bar: K));
    [Z_bar, Z_tilde] = deal(Z(1: K_bar, :), Z(1 + K_bar: K, :));
    [cost_tilde, ~, ~, ~] = fisher_cost(D_tilde', Y_range);
    cost_new = 0.5 * opts.lambda2 * cost_tilde;
    
    cost_diff = cost_old - cost_new;
end

function [cost_assigned] = cost_k(D, Y_range, K_bar, K_tilde, idx, opts)
    nClasses = 2;
    cost_assigned = [];
    y_k = D(:, idx)';
    % [~, p, ~, ~] = ttest2(y_k(1:Y_range(2))', y_k(1 + Y_range(2):Y_range(3))',.05);

    m_k = mean(y_k);
    % Assume D(:, idx) is assigned to D_bar
    bt = 0; wn = 0;
    for c = 1: nClasses
        y_k_c = get_block_col(y_k, c, Y_range);
        Nc = size(y_k_c, 2);
        m_ck = mean(y_k_c);
        bt = bt + Nc * (m_ck - m_k)^2;
        for n = 1: Nc
            wn = wn + (y_k_c(n) - m_ck)^2/Nc;
        end
    end
    for k = 1: K_bar
        cost_assigned = [cost_assigned, 0];
    end
    
    % Assume D(:, idx) is assigned to D_tilde
    cost = wn - bt + 2 * normF2(y_k);
    for k = 1 : K_tilde
        cost_assigned = [cost_assigned, opts.lambda2 * cost];       
    end
end
function [cost, Sw, Sb, cost_k] = fisher_cost(X, Y_range)
    nClasses = numel(Y_range) - 1;
    m = mean(X,2);
    n_atoms = size(X, 1);
    cost = 2 * normF2(X);
    Sw = 0;
    Sb = 0;
    for c = 1: nClasses
        Xc = get_block_col(X, c, Y_range);
        Mc = buildMean(Xc);
        Nc = size(Xc,2);
        M = repmat(m, 1, Nc);
        cost = cost + normF2(Xc - Mc) - normF2(Mc - M);
        Sw = Sw + normF2(Xc - Mc);
        Sb = Sb + normF2(Mc - M);
    end 
    cost_k = zeros(1, n_atoms);
    for k = 1: n_atoms
         for c = 1: nClasses
             Xc = get_block_col(X, c, Y_range);
             Nc = size(Xc,2);
             mc = mean(Xc, 2);
             for n = 1: Nc
                  cost_k(k) = cost_k(k) + (Xc(k, n) - mc(k))^2;
             end
             cost_k(k) = cost_k(k) - Nc * (mc(k) - m(k))^2;
         end
         cost_k(k) = cost_k(k) + normF2(X(k, :));
    end
end 