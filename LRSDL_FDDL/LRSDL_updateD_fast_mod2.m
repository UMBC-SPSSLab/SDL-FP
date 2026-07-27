function [D_bar, D_tilde, iter] = LRSDL_updateD_fast_mod2(Y, Y_range, D_bar, D_tilde, Z_bar, Z_tilde, opts)
    Z = [Z_bar; Z_tilde];
    A = Z*Z';
    B = Y*Z';
    [D_bar, D_tilde, iter] = ODL_updateD_mod2(D_bar, D_tilde, A, B, Y_range, opts);   
end