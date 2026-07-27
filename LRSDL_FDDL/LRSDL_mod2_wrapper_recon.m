function [] = LRSDL_mod2_wrapper_recon(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda1, lambda2, num_initial,seed,testindi,inittest)
% function [acc, rt] = LRSDL_mod2_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
% modified by S.-J.Kim for Asilomar'18

    if nargin == 0 % test mode
        %SK:
        error('nargin == 0 is not supported');
        
        dataset = 'myYaleB';
        N_train = 10;        
        [~, Y_train, Y_test, label_train, label_test] = ...
            train_test_split(dataset, N_train);        
        k = 8;
        k0 = 5;
        lambda1 = 0.001;
        lambda2 = 0.01;
        
    end 
    C = max(label_train);
    opts.k           = k;
    opts.k0          = k0;
    opts.show_cost   = 1;
    opts.lambda1     = lambda1;
    opts.lambda2     = lambda2;
    opts.D_range     = k*(0:1);
    opts.D_range_ext = [opts.D_range k+k0];  %SK:[opts.D_range k*C+k0];
    opts.initmode    = 'normal';  
    opts.max_iter    = 1000000;
    opts             = initOpts(opts);
    opts.verbose     = 0;
    opts.tol         = 1e-12;
    %% Train -- Note Y_train *transpose* is input here
    %%inittest=1, test different initializations, hence jump cost judgement
if inittest
    index_init=num_initial;
    [D, D0, X, X0, CoefM, coefM0, opts, rt, cost_final_old] = ...
                    LRSDL_mod2(Y_train', label_train, opts,index_init,seed,inittest);
     index_init_min=index_init;           
else
    for index_init=1:num_initial %initializers, different seeds have same local starting point.
    [D_new, D0_new, X_new, X0_new, CoefM_new, coefM0_new, opts_new, rt_new, cost_final_new] = ...
                    LRSDL_mod2(Y_train', label_train, opts,index_init,seed,inittest);
    %%Comparision of cost funciton values with different initializers and
    %%keep D,D0,X,X0,M,M0,opt,rt with minimum cost function value. 
    if index_init>=2
        if cost_final_old>cost_final_new
            D=D_new;
            D0=D0_new;
            X=X_new;
            X0=X0_new;
            CoefM=CoefM_new;
            coefM0=coefM0_new;
            opts=opts_new;
            rt=rt_new;
            cost_final_old=cost_final_new;
            index_init_min=index_init;
        end
    else
        D=D_new;
        D0=D0_new;
        X=X_new;
        X0=X0_new;
        CoefM=CoefM_new;
        coefM0=coefM0_new;
        opts=opts_new;
        rt=rt_new;
        cost_final_old=cost_final_new;
        index_init_min=index_init;%order of initializer for minimum value of cost function
    end  
    end
end

D_val = [X;X0]'\Y_test;
recon_err = norm(Y_test - [X;X0]'*D_val)^2;
sparsity = (sum(X(:) == 0) + sum(X0(:) == 0))/38/48546;
fn = fullfile("data", "recon_sparsity", num2str(lambda1) + "_" + num2str(lambda2) ...
    + "_" + num2str(sparsity) + "_" + num2str(recon_err) + "_.mat");
save(fn, "lambda1", "lambda2", "sparsity", "recon_err");
