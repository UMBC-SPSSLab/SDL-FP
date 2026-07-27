function [Y_train, label_train, Y_val, label_val] = pickTrainTest_3(dataset, N_train_c, N_val_c, N_test_c, seed)
    rng(seed, 'twister');
    data_fn = fullfile('data', strcat(dataset, '.mat'));
    load(data_fn);

    Y = normc(Y);
    d = size(Y,1);
    if ~exist('Y_range', 'var')
        Y_range = label_to_range(label);
    end

    C = numel(Y_range) - 1;
    N_train = C*N_train_c;
    N_val=C*N_val_c;
    N_test=C*N_test_c;
    Y_train = zeros(d, N_train);
    Y_val = zeros(d, N_val);
    Y_test=zeros(d, N_test);
    label_train = zeros(1, N_train);
    label_val = zeros(1, N_val);
    label_test = zeros(1, N_test);

    cur_train = 0;
    cur_val = 0;
    cur_test = 0;
    for c = 1: C 
        Yc = get_block_col(Y, c, Y_range);
        N_total_c = size(Yc, 2);
        label_train(:, cur_train + 1: cur_train + N_train_c) = c*ones(1, N_train_c);
        label_val(:, cur_val + 1: cur_val + N_val_c) = c*ones(1, N_val_c);
        label_test(:, cur_test + 1: cur_test + N_test_c) = c*ones(1, N_test_c);
        idx = randperm(N_total_c);

        Y_train(:, cur_train + 1: cur_train + N_train_c) = Yc(:, idx(1: N_train_c));
        Y_val(:, cur_val + 1: cur_val + N_val_c) = Yc(:, idx(N_train_c + 1: N_train_c +N_val_c));
        Y_test(:, cur_test + 1: cur_test + N_test_c) = Yc(:, idx(N_train_c + N_val_c+1: N_train_c +N_val_c+N_test_c));

        cur_train = cur_train + N_train_c;
        cur_val = cur_val + N_val_c;
        cur_test= cur_test+N_test_c;
    end 


    Y_train = normc(Y_train);
    Y_val = normc(Y_val);
    Y_test = normc(Y_test);
end 


