function [acc_knn, acc_ksvm, D_val] = LRSDL_mod2_pred_GC(Y_val, Z_bar, Z_tilde, D_train_tilde, opts, label_val, label_train)
    k_tilde = opts.k_tilde;
    k_bar = opts.k_bar;

    D_range = k_tilde*(0:1);
    
    N = size(Y_val, 2);
    
    D_val = [Z_bar; Z_tilde]'\Y_val;
    D_val = D_val';
    D_val_tilde = D_val(:, k_bar + 1: end);
    %{
    for k = 1: k_tilde
         D_val_tilde(:, k) = D_val_tilde(:, k) / max(1, norm(D_val_tilde(:, k)));
    end
    
    %% --------- classification by nearest centroid -------
    nClasses = 2;
    D_train_tilde_Mean = zeros(size(D_train_tilde, 2), nClasses);
    subjects_range_train = label_to_range(label_train);
    for c = 1: nClasses
        Xc = get_block_col(D_train_tilde', c, subjects_range_train);
        D_train_tilde_Mean(:, c) = mean(Xc, 2); % mean of each row(atoms)
    end
    E = zeros(nClasses, N);
    for c = 1: nClasses            
        Mc = repmat(D_train_tilde_Mean(:, c), 1, N);
        R1 = (D_val_tilde' - Mc);
        E(c,:) = sum(R1.^2);
    end
    [~, pred] = min(E);
    acc_NC = double(sum(pred == label_val))/N;
    fprintf('opts.lambda1 = %.4f,  acc_NC: %f\n', opts.lambda1, acc_NC);
    %}
    %% ----- alternative classificaiton by knn ----------------
    Mdl = ExhaustiveSearcher(D_train_tilde);
    [idx,~] = knnsearch(Mdl, D_val_tilde, 'K', 1);
    pred = label_train(idx);
    acc_knn = double(sum(pred == label_val))/N;
    
    %{%}
    svmModel = fitcsvm(D_train_tilde, label_train, 'KernelFunction', 'linear', 'BoxConstraint', 1);
    predictedLabels = predict(svmModel, D_val_tilde)';
    acc_ksvm = double(sum(predictedLabels == label_val))/N;
    
      
    
    fprintf('opts.lambda1 = %.4f,  acc_knn: %f, acc_ksvm: %f\n', opts.lambda1, acc_knn, acc_ksvm);
end 