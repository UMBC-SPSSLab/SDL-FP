 function training()
    addpath('utils');
    addpath('LRSDL_FDDL');
    addpath('ODL');
    dataset = "AOD1";
    isTuning = 0;
    range = 75;
    n = 75;
    seed = 4;
    rng(seed,'twister');
    random_ints = [1:100]; % sort(randperm(range, n));
    % read parameters
    fid = fopen("data/bestPred_parameters.txt", "r"); 
    format = "%f %f %f %d %d %d %d %d %d";
    size_pars = [9 Inf];
    pars = fscanf(fid, format, size_pars);
    [lambda1, lambda2, lambda3, K_bar, K_tilde, isJSTSP, isSimulation, permutate, total_run] = deal(pars(1, 1), pars(2, 1), pars(3, 1), pars(4, 1), pars(5, 1), pars(6, 1), pars(7, 1), pars(8, 1), pars(9, 1));
    fclose(fid);
    if 1 == isSimulation
        dataset = "Synthesized_dataset";
    end
    % load data and Y is voxels by subjects
    
    data_fn = fullfile('data', strcat(dataset, '.mat'));
    load(data_fn);
    Y = normc(Y);
    tic;
    while true
        %%%%%%%%%%%%%%%%%%%% read next task ID %%%%%%%%%%%%%%%%%%
        fid = fopen("data/validRuns.txt", "r+");
        format = "%d";
        read_data = fscanf(fid, format);
        seedID = read_data(end); 
        if seedID > total_run
            break;
        end 
        fprintf(fid, "\n%d", seedID + 1);
        fclose(fid);
        seedID = random_ints(seedID);
        
        fprintf('Starting. lambda1 = %f, lambda2 = %f, lambda3 = %f, K_bar = %d, K_tilde = %d, seedID = %d\n', lambda1, lambda2, lambda3, K_bar, K_tilde, seedID);
        t1 = toc();                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 label = range_to_label(Y_range);
        [D_bar, D_tilde, Z_bar, Z_tilde, ~, ~, ~, recon, cost, costs] = LRSDL_solveDZ(Y, label, [], [], K_tilde, K_bar, lambda1, lambda2, lambda3, seedID, isTuning, isJSTSP, permutate);
        t2 = toc;
        fprintf("Training took %f seconds\n", t2 - t1);
        mkdir("data", "results" + "/" + string(seedID));
        fn = fullfile("data", "results", string(seedID), "D_bar.mat");
        save(fn, "D_bar");    
        fn = fullfile("data", "results", string(seedID), "D_tilde.mat");
        save(fn, "D_tilde");    
        fn = fullfile("data", "results", string(seedID), "Z_bar.mat");
        save(fn, "Z_bar");   
        fn = fullfile("data", "results", string(seedID), "Z_tilde.mat");
        save(fn, "Z_tilde");   
        fn = fullfile("data", "results", string(seedID), "recon.mat");
        save(fn, "recon");  
        fn = fullfile("data", "results", string(seedID), "info.mat");
        save(fn, "lambda1", "lambda2", "lambda3", "Y_range", "costs", "cost", "K_tilde", "K_bar");
    end
end