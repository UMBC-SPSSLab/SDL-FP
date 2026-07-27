function [D_bar, D_tilde] = LRSDL_updateD_fast_mod2_lr(Y, Y_range, D_bar, D_tilde, Z_bar, Z_tilde, opts)    
    Y_tilde = Y - D_bar * Z_bar;
    D_tilde = ODL_updateD_fisher_lr(D_tilde, Z_tilde * Z_tilde', Y_tilde * Z_tilde', Y_range, opts);
    
    Y_bar = Y - D_tilde * Z_tilde;
    D_bar = min_rank_dict(D_bar, Y_bar * Z_bar', Z_bar * Z_bar', opts.lambda3, opts);   
end