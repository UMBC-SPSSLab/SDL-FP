function tuning()
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    [N_train, N_val, N_test] = deal(85, 36, 0);
    isTuning = 1;
    % dataset = "AOD1";
    dataset = "Synthesized_dataset";

    numInit = 3;
    numFold = 3;
    tic;
    while true
        % find how many combinations of parameters are there
        fid = fopen("data/tasks_parameters.txt", "r"); 
        format = "%d %f %f %f %d %d %d %d %d";
        size_pars = [9 Inf]; 

        pars = fscanf(fid, format, size_pars);
        num_par = size(pars, 2);
        fclose(fid);
        %%%%%%%%%%%%%%%%%%%% read next task ID %%%%%%%%%%%%%%%%%%
        fid = fopen("data/tasks.txt", "r+");
        format = "%d";
        parIDs = fscanf(fid, format);
        parID = parIDs(end); 
        if parID > num_par
            break;
        end 
        fprintf(fid, "\n%d", parID + 1);
        fclose(fid);
        [init_seed, lambda1, lambda2, lambda3, k_bar, k_tilde, isJSTSP, isSimulation, permutate] = ...
                deal(pars(1, parID), pars(2, parID), pars(3, parID), pars(4, parID), pars(5, parID), pars(6, parID), pars(7, parID), pars(8, parID), pars(9, parID));
        if 1 == isSimulation
            dataset = "Synthesized_dataset";
        end
        fprintf("seed = %d, lambda1 = %f, lambda2 = %f, lambda3 = %f, k_bar = %d, k_tilde = %d, isJSTSP = %d \n", init_seed, lambda1, lambda2, lambda3, k_bar, k_tilde, isJSTSP);
        
        fprintf('starting... isJSTSP %d, %s\n', isJSTSP, dataset);
        spmd(numFold)
            [accs_knn, accs_ksvm, recon_errors, costs] = deal(0, 0, 0, 0);
            % Y_train/val are voxels by subjects
            [Y_train, label_train, Y_val, label_val] = pickTrainTest_2(dataset, N_train, N_val, N_test, init_seed + labindex);
            for i = 1 : numInit
                [D_bar, D_tilde, Z_bar, Z_tilde, D_val, acc_knn, acc_ksvm, recon_error, cost, ~] = LRSDL_solveDZ(Y_train, label_train, Y_val, label_val, ...
                                k_tilde, k_bar, lambda1, lambda2, lambda3, init_seed + i, isTuning, isJSTSP, permutate);
                [accs_knn, accs_ksvm, recon_errors, costs] = deal(accs_knn + acc_knn, accs_ksvm + acc_ksvm, recon_errors + recon_error, costs + cost);
                %mySave(i, labindex, D_bar, D_tilde, Z_bar, Z_tilde, label_train, label_val, acc_knn, acc_ksvm, D_val);
            end
        end
        [acc_knn, acc_ksvm, recon_error, cost] = deal(0, 0, 0, 0);
        for i = 1: numFold
            [acc_knn, acc_ksvm, recon_error, cost] = deal(acc_knn + accs_knn{i}, acc_ksvm + accs_ksvm{i}, recon_error + recon_errors{i}, cost + costs{i});
        end
        [acc_knn, acc_ksvm, recon_error, cost] = deal(acc_knn/numInit/numFold, acc_ksvm/numInit/numFold, recon_error/numInit/numFold, cost/numInit/numFold);
        fprintf("Done ! seed = %d, lambda1 = %f, lambda2 = %f, lambda3 = %f, k_bar = %d, k_tilde = %d, isJSTSP = %d \n", init_seed, lambda1, lambda2, lambda3, k_bar, k_tilde, isJSTSP);
        fn = fullfile("data", "results", string(lambda1) + "_" + string(lambda2) + "_" + string(lambda3) + "_" + ...`
            string(acc_knn) + "_" + string(acc_ksvm) + "_" + string(recon_error) + "_" + string(cost) + "_" + string(init_seed) + "_.mat");
        save(fn, "lambda1");    
    end
    disp("Tunning done, took:");
    disp(toc/3600);
    disp("Hours");
end
function mySave(id_init, labindex, D_bar, D_tilde, Z_bar, Z_tilde, label_train, label_val, acc_knn, acc_ksvm, D_val)
    mkdir("data", "results" + "/" + string(labindex) + "_" + string(id_init));
    fn = fullfile("data", "results" + "/" + string(labindex) + "_" + string(id_init), "results.mat");
    save(fn, "D_bar", "D_tilde", "Z_bar", "Z_tilde", "label_train", "acc_knn", "acc_ksvm", "D_val", "label_val");
end

