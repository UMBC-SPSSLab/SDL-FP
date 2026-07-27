filePath = "..\data\ResultsofDLandICA";
fileNames = ["DLbestrun_Common_Comp.mat", "DLbestrun_Disc_Comp.mat", "DLbestrun_Common_Act.mat", "DLbestrun_Disc_Act.mat"];

K_bar_J = 12;
K_tilde_J = 26;

for i = 1: size(fileNames, 2)
    fn = fullfile(filePath, fileNames(i)); 
    load(fn);
end
Z_J = [X0; X]; % maps (components) of JSTSP
D_bar_J = D0;
D_tilde_J = D;
Y_range_J = [0, 121, 242];
clear D0 D;
filePath = "..\experiments\2NormY\real-permutate\bar_21_tilde_17\0.003 0.52 0.55 - NoPert\35";
fileNames = ["Z_bar.mat", "Z_tilde.mat", "info.mat", "D_bar", "D_tilde"];
for i = 1: size(fileNames, 2)
    fn = fullfile(filePath, fileNames(i)); 
    load(fn);
end
K_bar_P = 21;
K_tilde_P = 17;
D_bar_P = D_bar;
D_tilde_P = D_tilde;
Y_range_P = Y_range;
clear D_bar D_tilde Y_range;

Z_P = [Z_bar; Z_tilde]; % maps (components) of proposed method

K = 38;
p2j_corr = zeros(K, K);
for i = 1:K
    for j = 1:K
        tmp = corrcoef(Z_P(i, :)', Z_J(j, :)');
        p2j_corr(i, j) = tmp(1,2);
    end
end 
[mx_cor, idxs] = max(abs(p2j_corr)');
cor = [mx_cor; idxs];
cor_bar = cor(:, 1 : K_bar_P);
cor_tilde = cor(:, 1 + K_bar_P: K);

for k = 1: K
    mx_cor(k) = p2j_corr(k, idxs(k)); 
end
figure; imagesc((p2j_corr));


% [proposed, JSTSP]
p2j_corr = -1 * abs(p2j_corr);
assign_proposed_JSTSP = matchpairs(p2j_corr, 1000000);

j2p_corr = p2j_corr';
assign_JSTSP_proposed = matchpairs(j2p_corr, 1000000);

% check if there are two same maps in proposed results
self_corr_proposed = zeros(K, K);
for i = 1:K
    for j = 1:K
        tmp = corrcoef(Z_P(i, :)', Z_P(j, :)');
        self_corr_proposed(i, j) = tmp(1,2);
    end
end
clear tmp;
figure; imagesc(self_corr_proposed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Now Plot The Matched Maps %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('indices.mat');
addpath('spm12'); %spm:statistical parametric mapping
anat = spm_read_vols(spm_vol('nsingle_subj_T1_2_2_3.hdr'));
z_choice = 4:5:46;

Sf_J = Z_J;
Sf_P = Z_P;

dat3f_J = zeros(K,53*63*46); 
dat3f_J(:,index_) = Sf_J(:,(1:length(index_)));
dat3f_J = reshape(dat3f_J, K, 53, 63, 46);

dat3f_P = zeros(K,53*63*46); 
dat3f_P(:,index_) = Sf_P(:,(1:length(index_)));
dat3f_P = reshape(dat3f_P, K, 53, 63, 46);

%%%%%%%%%% Get the sign ambiguity %%%%%%

p_bar_P = zeros(1, K_bar_P);
p_tilde_P = zeros(1, K_tilde_P);
p_bar_J = zeros(1, K_bar_J);
p_tilde_J = zeros(1, K_tilde_J);
T_P = zeros(1, K);
T_J = zeros(1, K);
sele_P = zeros(1, K);
sele_J = zeros(1, K);
monthresh = 2.5; %set z-score threshold (higher threshold means fewer plotted voxels

idx = 1;
for k = 1: K_bar_P
    [~, p_bar_P(idx), ~, tstats] = ttest2(D_bar_P(1:Y_range_P(2),k), D_bar_P(1 + Y_range_P(2):Y_range_P(3),k),.05);
    T_P(idx)=tstats.tstat;
    sele_P(idx) = sign(T_P(idx))*idx;
    idx = idx + 1;
end

for k = 1: K_tilde_P
    [~, p_tilde_P(idx - K_bar_P), ~, tstats] = ttest2(D_tilde_P(1:Y_range_P(2),k), D_tilde_P(1 + Y_range_P(2):Y_range_P(3),k),.05);
    T_P(idx) = tstats.tstat;
    sele_P(idx) = sign(T_P(idx))*idx;
    idx = idx + 1;
end

idx = 1;
for k = 1: K_bar_J
    [~, p_bar_J(idx), ~, tstats] = ttest2(D_bar_J(1:Y_range_J(2),k), D_bar_J(1 + Y_range_J(2):Y_range_J(3),k),.05);
    T_J(idx)=tstats.tstat;
    sele_J(idx) = sign(T_J(idx))*idx;
    idx = idx + 1;
end

for k = 1: K_tilde_J
    [~, p_tilde_J(idx - K_bar_J), ~, tstats] = ttest2(D_tilde_J(1:Y_range_J(2),k), D_tilde_J(1 + Y_range_J(2):Y_range_J(3),k),.05);
    T_J(idx) = tstats.tstat;
    sele_J(idx) = sign(T_J(idx))*idx;
    idx = idx + 1;
end

%{
% 4 figures, up to 10 maps each figure
colums = 10;
for num_fig = 1: 4
    figure;
    for num_col = 1: colums
        idx = (num_fig - 1) * colums + num_col;
        if idx > K
            break; 
        end
        % first row: proposed map
        tmp(:, :, :) = dat3f_P(abs(sele_P(idx)),:,:,z_choice);%x,y,time,slice?
        subplot_tight(2, colums, num_col, [0.001,0.005])
        xlabel('fMRI', 'FontSize', 12);
        colorbar;
        make_composite(sign(sele_P(idx))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component

        % second row: JSTSP map
        idx = assign_JSTSP_proposed(idx, 1); % index of JSTSP map
        tmp(:, :, :) = dat3f_J(abs(sele_J(idx)),:,:,z_choice);%x,y,time,slice?
        subplot_tight(2, colums, num_col + colums, [0.001,0.005])
        xlabel('fMRI', 'FontSize', 12);
        colorbar;
        make_composite(sign(sele_J(idx))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component
    end
end
%}
% indexes of maps in proposed results that are improved
improved_P = [];
% indexes of maps in proposed results that are consistent with JSTSP
% results
same_P = [];
different_P = [];
unmatched_P = [];
p_j = [p_bar_J, p_tilde_J];
p_p = [p_bar_P, p_tilde_P];

for i_p = 1: K_bar_P
    i_j = assign_JSTSP_proposed(i_p, 1);
    if abs(j2p_corr(i_j, i_p)) > 0.3
        if i_j <= K_bar_J
            if p_j(i_j) >= 0.05
                same_P = [same_P, i_p];
            else
                improved_P = [improved_P, i_p];
            end
        else
            if p_j(i_j) >= 0.05
                improved_P = [improved_P, i_p];
            else
                different_P = [different_P, i_p];
            end
        end
    else
        unmatched_P = [unmatched_P, i_p];
    end
end

for i_p = 1 + K_bar_P: K_bar_P + K_tilde_P
    i_j = assign_JSTSP_proposed(i_p);
    if abs(j2p_corr(i_j, i_p)) > 0.3
        if i_j <= K_bar_J
            if p_j(i_j) >= 0.05
                different_P = [different_P, i_p];
            else
                improved_P = [improved_P, i_p];
            end
        else
            if p_j(i_j) >= 0.05
                improved_P = [improved_P, i_p];
            else
                same_P = [same_P, i_p]; 
            end
        end
    else
        unmatched_P = [unmatched_P, i_p];
    end
end

new_permation = [improved_P, same_P, different_P, unmatched_P];
new_permation_j = assign_JSTSP_proposed(new_permation, 1);

% 7 figures, up to 9 maps each figure
colums = 9;
for num_fig = 1: 5
    figure;
    for num_col = 1: colums
        idx = (num_fig - 1) * colums + num_col;
        if idx > K
            break; 
        end
        % first row: proposed map
        tmp(:, :, :) = dat3f_P(abs(sele_P(new_permation(idx))),:,:,z_choice);%x,y,time,slice?
        subplot_tight(2, colums, num_col, [0.001,0.005])
        xlabel('fMRI', 'FontSize', 12);
        colorbar;
        make_composite(sign(sele_P(new_permation(idx)))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component

        % second row: JSTSP map
        idx = assign_JSTSP_proposed(new_permation(idx), 1); % index of JSTSP map
        tmp(:, :, :) = dat3f_J(abs(sele_J(idx)),:,:,z_choice);%x,y,time,slice?
        subplot_tight(2, colums, num_col + colums, [0.001,0.005])
        xlabel('fMRI', 'FontSize', 12);
        colorbar;
        make_composite(sign(sele_J(idx))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component
    end
end
matched_corr_p2j = zeros(1, K);
for i = 1: K
    matched_corr_p2j(i) = p2j_corr(i, assign_JSTSP_proposed(i, 1));
end

for i = 1: 38
    fprintf("%.4f           ", p_p(new_permation(i)));
    if i == 9 || i == 18 || i == 27 || i == 36
        fprintf("\n");
    end
end
fprintf("\n---------\n");
for i = 1: 38
    fprintf("%.4f           ", p_j(assign_JSTSP_proposed(new_permation(i), 1)));
    if i == 9 || i == 18 || i == 27 || i == 36
        fprintf("\n");
    end
end
clear; close all;