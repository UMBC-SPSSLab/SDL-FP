% Top level code for choosing best DL run from all the given runs
addpath_recurse("GroupICATv4.0c");
path = "experiments\2NormY\simulation-permutate\step1.5\datasets_10001\jstsp_per\results\runs";
% path = "data\results";
files = dir(path);
files = files(3:end, :);
numericNames = str2double({files.name});
[~, sortIdx] = sort(numericNames);
files = files(sortIdx);
no_runs = length(files);
run_idxes = [];
Zs = {};

for i = 1: length(files)
    run_idx = files(i).name;
    run_idxes = [run_idxes, str2num(run_idx)];
    data_fn = fullfile(path, run_idx, 'Z_tilde.mat');
    load(data_fn);
    data_fn = fullfile(path, run_idx, 'Z_bar.mat');
    load(data_fn);
    Zs{i} = [Z_bar; Z_tilde];
end 
[StdT, RunCorrT, BestZ, BestZIndex, IndexAssign, Tmaps] = DL_bestRunSelection(Zs);

fprintf("BestRunIndex = %d\n", run_idxes(BestZIndex));

sum_RunCorrT = sum(RunCorrT, 2);
[~,idx] = sort(sum_RunCorrT);
RunCorrT = RunCorrT(flip(idx), :);

figure; imagesc(RunCorrT); colorbar; title("RunCorrT Proposed");

figure; 
cdfplot(StdT(:)); title("CDF StdT Proposed");

clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% For JSTSP %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
path = "data\ResultsofDLandICA\N_121_K_26_K0_12_l1_0.003_l2_0.2";
files = dir(path);
Zs = {};
costs = [];
recons = [];
[s, e] = deal(1, 100);
for i = s + 2: length(files) - 1
    sp = strsplit(files(i).name, '_');
    z = sp{12};
    a = sp{13};
    b = sp{17};
    file_path = strcat('data\ResultsofDLandICA\N_121_K_26_K0_12_l1_0.003_l2_0.2\AOD1_N_121_k_26_k0_12_l1_0.003_l2_0.2_', z, '_', a, '_seed_50_initorder_', b);
    data_fn = fullfile(file_path);
    load(data_fn);
    sp = strsplit(b, '.');
    num = str2num(sp{1});
    if num <= e
        Zs{num} = [X; X0];
    end
end 
[StdT_JSTSP, RunCorrT_JSTSP, BestZ, BestZIndex_JSTSP, IndexAssign, Tmaps_JSTSP] = DL_bestRunSelection(Zs);

%%%%%%%%% compute correlations of runs w.r.t. T-maps %%%%%%%%%
sum_RunCorrT_JSTSP = sum(RunCorrT_JSTSP, 2);

%%% sort each column in increasing order %%%%
[out,idx] = sort(sum_RunCorrT_JSTSP);

%%% rearrange rows resulting in each column in decreasing order %%%%
RunCorrT_JSTSP = RunCorrT_JSTSP(flip(idx), :);

figure; imagesc(RunCorrT_JSTSP); colorbar; title("RunCorrT JSTSP");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[max_Tmaps_JSTSP, ~] = max(abs(Tmaps_JSTSP)'); 
cdfplot(max_Tmaps_JSTSP(:)); 
hold on;
[max_Tmaps, ~] = max(abs(Tmaps)'); 
cdfplot(max_Tmaps(:)); 
title("max\_Tmaps");
legend('max\_Tmaps\_JSTSP', 'max\_Tmaps\_Proposed', 'best', 'Location', "east");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
[max_StdT_JSTSP, ~] = max(abs(StdT_JSTSP)'); 
cdfplot(max_StdT_JSTSP(:)); 
hold on;
[max_StdT, ~] = max(abs(StdT)'); 
cdfplot(max_StdT(:)); 
title("max\_StdT");
legend('JSTSP', 'Proposed', 'best', 'Location', "east");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
mean_StdT_JSTSP = mean(abs(StdT_JSTSP)'); 
cdfplot(mean_StdT_JSTSP(:)); 
hold on;
mean_StdT = mean(abs(StdT)'); 
cdfplot(mean_StdT(:)); 
title("mean\_StdT");
legend('JSTSP', 'Proposed', 'best', 'Location', "east");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
cdfplot(abs(StdT_JSTSP(:))); 
hold on;
cdfplot(abs(StdT(:))); 
title("cdf\_StdT");
legend('JSTSP', 'Proposed', 'best', 'Location', "east");
clear; close all;
%}