function [acc,D_test,D0_test] = LRSDL_mod2_test_wrapper(Y_test,X,X0,CoefM,coefM0,D,D0, ...
                            k, k0, lambda1, lambda2)
% function [acc, rt] = LRSDL_mod2_wrapper(Y_train, label_train, Y_test , label_test, ...
%                             k, k0, lambda1, lambda2)
% -----------------------------------------------
% Author: Tiep Vu, thv102@psu.edu, 5/11/2016
%         (http://www.personal.psu.edu/thv102/)
% -----------------------------------------------
% modified by S.-J.Kim for Asilomar'18
label_train=[ones(1,60) ones(1,60)*2];
label_test=[ones(1,30) ones(1,30)*2];

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
    opts.max_iter    = 100;
    opts             = initOpts(opts);
    opts.verbose     = true;
    opts.tol         = 1e-8;
    
    %% Train -- Note Y_train *transpose* is input here
%     [D, D0, X, X0, CoefM, coefM0, opts, rt] = ...
%                     LRSDL_mod2(Y_train', label_train, opts);
%     
    %% Testing
%     X1 = [X; X0];
%     Y_range = label_to_range(label_train);
%     C = max(label_train);
%     CoefMM0 = zeros(size(X1,1), C);
%     for c = 1: C 
%         X1c = get_block_col(X1, c, Y_range);
%         CoefMM0(:,c) = mean(X1c,2);
%     end    
    opts.verbose = 0;
    acc = [];
    if numel(D0) ~= 0
        fprintf('GC:\n');
        %SK: acc1 = LRSDL_pred_GC(Y_test, D, D0, CoefM, coefM0, opts, label_test);
        [acc1,D_test,D0_test] = LRSDL_mod2_pred_GC(Y_test, X, X0, CoefM, coefM0, opts, label_test, ...
            D, label_train);
%         fprintf('LC:\n');
%         acc2 = LRSDL_pred_LC(Y_test, D, D0, CoefMM0, opts, label_test);
%         acc = [acc1 acc2];
        fprintf('maximum acc: %4f\n', max(acc));
    else             
        error('not supported');
        fprintf('GC:\n');
        opts.weight = 0.5;
        for vgamma = [0.0001, 0.001, 0.01, 0.1]
            opts.gamma = vgamma;
            pred = FDDL_pred(Y_test, D, CoefM, opts);
            acc = [acc calc_acc(pred, label_test)];
            fprintf('gamma = %.4f, acc = %.4f\n', vgamma, acc(end));
        end 
        %% LC Uncomment following line to perform LC prediction.
%         fprintf('LC:\n');
%         gamma2 = 0.01;
%         opts.gamma2 = gamma2;
%         for gamma1 = [0.0001, 0.001, 0.01, 0.1]
%             opts.gamma1 = gamma1;
%             pred = FDDL_pred_LC(Y_test, D, CoefM, opts);
%             acc = [acc calc_acc(pred, label_test)];
%             fprintf('gamma = %.4f, acc = %.4f\n', gamma1, acc(end));
%         end 
    end
    fprintf('\n');
    acc = max(acc1);
end 
