function [D_bar, D_tilde, Z_bar, Z_tilde, D_val, acc, error_fit, cost_final_old] = LRSDL_solveDZ_lr(Y_train, label_train, Y_val , label_val, ...
                            k_tilde, k_bar, lambda1, lambda2, lambda3, init_seed, isTuning, isJSTSP)
    C = max(label_train);
    opts.k_tilde    = k_tilde;
    opts.k_bar      = k_bar;
    opts.lambda1    = lambda1;
    opts.lambda2    = lambda2;
    opts.lambda3    = lambda3;
    opts.isJSTSP = isJSTSP;
    opts.D_col_range = k_tilde*(0:1);
    opts.D_col_range_ext = [opts.D_col_range k_tilde + k_bar]; 
    opts.initmode   = 'normal';  
    opts.max_iter_D   = 2000;
    opts.max_iter   = 300;
    opts            = checkOpts(opts);
    opts.verbose    = 0;
    opts.verboseD    = 0;
    opts.tol        = 1e-12;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Training, find a decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [D_bar, D_tilde, Z_bar, Z_tilde, opts, cost_final_old] = LRSDL_mod2_lr(Y_train', label_train, opts, init_seed);
                
    %%%%%%%%%%%%%%%%%%%%%%%%% Tuning, compute prediction accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%
    [error_fit, acc, D_val] = deal(0, 0, []);
    if isTuning
        [acc, D_val] = LRSDL_mod2_pred_GC(Y_val, Z_bar, Z_tilde, D_tilde, opts, label_val, label_train);
        Z = [Z_bar; Z_tilde];
        error_fit = normF2(Y_val' - D_val*Z);
    end
end
