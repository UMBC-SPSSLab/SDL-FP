 function [cost_rec, cost_plain, cost_fisher, cost_inverFisher, cost] = LRSDL_cost_mod2(Y, Y_range, D_bar, D_tilde, Z_bar, Z_tilde, opts)
    nClasses = numel(Y_range) - 1;
 
    cost_rec = .5*normF2(Y - [D_bar, D_tilde] * [Z_bar; Z_tilde]); 
    cost_plain = opts.lambda1 * norm1([Z_bar; Z_tilde]);
    cost_fisher = 0.5 * opts.lambda2 * FDDL_discriminative(D_tilde', Y_range);

    cost_inverFisher = 0;
    if ~opts.isJSTSP
        cost_inverFisher = 0.5 * opts.lambda3 * trace(D_bar' * opts.H_bar * D_bar);
    end
    cost = cost_rec + cost_plain + cost_fisher + cost_inverFisher;
end 