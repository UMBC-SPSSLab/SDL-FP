% filePath = "..\data\results\2";

filePath = "..\experiments\2NormY\dict_bilevel_opt\real\v6_3\22 16";
fn = fullfile(filePath, 'tuning_result_seed1.mat'); 
load(fn);

fn = fullfile(strcat(filePath, '/runs_seed1'), 'run_58.mat'); 
load(fn);

% D_bar = D0; D_tilde = D; Z_bar = X0; Z_tilde = X;
% Y_range = [0, 121, 242];  
K_bar = K_com;
K_tilde = K_dis;
K = K_bar + K_tilde;

D_bar = D_trained(:, 1:K_bar);
D_tilde = D_trained(:, K_bar + 1:end);
Z_bar = Z_trained(1:22, :);
Z_tilde = Z_trained(23:end, :);
Y_range = [0, 150, 271];


M = size(D_bar, 1);
V = size(D_bar, 2);

% Y_range = [0, 150, 271];

Z = [Z_bar; Z_tilde];
D = [D_bar, D_tilde];
Sf = Z;
load('indices.mat');
addpath('spm12');%spm:statistical parametric mapping
fileName=['test']; %filename for saving the figures
monthresh = 2.5; %set z-score threshold (higher threshold means fewer plotted voxels
sele = [1: K];
 
%% Create 4-D fMRI data and read in anatomical map
dat3f = zeros(K,53*63*46); 
dat3f(:,index_) = Sf(:,(1:length(index_)));
dat3f = reshape(dat3f, K, 53, 63, 46);% why use reshape here?

fprintf('- Reading in Anatomical Data\n');
anat = spm_read_vols(spm_vol('nsingle_subj_T1_2_2_3.hdr'));
%anat is the graph or image of anatomical sturcture?
% select slices to show ("1:1:46" will show all 46 slices)
z_choice = 4:5:46;

figure;
imagesc(D_bar);
title("Imagesc of D\_bar");
figure;
imagesc(D_tilde);
title("Imagesc of D\_tilde");

% get the T values from two-sample t-test
p_bar = zeros(1, K_bar);
p_tilde = zeros(1, K_tilde);
p_norm_bar = zeros(1, K_bar);
p_norm_tilde = zeros(1, K_tilde);
idx = 1;
for k = 1: K_bar
    [h(idx), p_bar(idx), ~, tstats(k)] = ttest2(D_bar(1:Y_range(2),k), D_bar(1 + Y_range(2):Y_range(3),k),.05);
    p_norm_bar(k) = norm(D_bar(:, k));
    T(idx)=tstats(k).tstat;
    sele(idx) = sign(T(idx))*idx;
    idx = idx + 1;
end

for k = 1: K_tilde
    [h(idx), p_tilde(idx - K_bar), ~, tstats(k)] = ttest2(D_tilde(1:Y_range(2),k), D_tilde(1 + Y_range(2):Y_range(3),k),.05);
    T(idx) = tstats(k).tstat;
    p_norm_tilde(k) = norm(D_tilde(:, k));
    sele(idx) = sign(T(idx))*idx;
    idx = idx + 1;
end

% plot Z_bar
pp_bar = -100 + zeros(ceil(K_bar/10), 10);
pp_norm_bar = -100 + zeros(ceil(K_bar/10), 10);

figure;
for s = 1: ceil(K_bar/10)
    for j = 1: 10
        i = (s - 1) * 10 + j;
        if i <= K_bar
            pp_bar(s, j) = p_bar(i);
            pp_norm_bar(s, j) = p_norm_bar(i);
            tmp(:, :, :) = dat3f(abs(sele(i)),:,:,z_choice);%x,y,time,slice?
            subplot_tight(ceil(K_bar/10), 10, i, [0.001,0.005]);
            xlabel('fMRI', 'FontSize', 12);
            colorbar;
            make_composite(double(sign(sele(i)))*tmp/stdN(tmp),anat(:,:,z_choice),monthresh); %plot component
            % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
            clear tmp
        end
    end
end 
 
% plot Z_tilde
pp_tilde = -100 + zeros(ceil(K_tilde/10), 10);
pp_norm_tilde = -100 + zeros(ceil(K_tilde/10), 10);

figure;
for s = 1: ceil(K_tilde/10)
    for j = 1: 10
        i = (s - 1) * 10 + j;
        if i <= K_tilde
            pp_tilde(s, j) = p_tilde(i);
            pp_norm_tilde(s, j) = p_norm_tilde(i);
            tmp(:,:,:) = dat3f(abs(sele(K_bar + i)), :, :, z_choice);%x,y,time,slice?
            subplot_tight(ceil(K_tilde/10), 10, i, [0.001,0.005]);
            xlabel('fMRI', 'FontSize', 12);
            colorbar
            make_composite(double(sign(sele(K_bar + i)))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component
            % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
            clear tmp
        end
    end
end 

incoher_cross = normF2(D_bar'*D_tilde)/K_bar/K_tilde;
incoher_bar = normF2(D_bar' * D_bar)/K_bar/K_bar;
incoher_tilde = normF2(D_tilde' * D_tilde)/K_tilde/K_tilde;
tmp = D'*D;
figure; imagesc(tmp); colorbar;
fprintf("incoher_cross = %f, incoher_bar = %f, incoher_tilde = %f\n", incoher_cross, incoher_bar, incoher_tilde);
xlabel("Incoherence");
figure;
cdfplot(abs(Z_bar(:)));
title("cdf of Z\_bar");

figure;
cdfplot(abs(Z_tilde(:)));
title("cdf of Z\_tilde");

geomean_bar = geomean(p_bar, 'all');
geomean_tilde = geomean(p_tilde, 'all');
geomedian_bar = geometric_median(p_bar);
geomedian_tilde = geometric_median(p_tilde);

fprintf("Geometric mean of p values of D_bar: %f, D_tilde: %f.\n", geomean_bar, geomean_tilde);
fprintf("Geometric median of p values of D_bar: %f, D_tilde: %f.\n", geomedian_bar, geomedian_tilde);
disp("pp_bar");
disp(pp_bar);
disp("--------------------");
disp("pp_tilde");
disp(pp_tilde);

pp = [pp_bar; [1:10] ;pp_tilde];
% [D, permutation] = solve_permutation(D_bar, D_tilde, K_bar, K_tilde, Y_range);

D_HC = D(1:Y_range(2), :);
D_SZ = D(1 + Y_range(2):Y_range(3), :);
mean_diff = mean(D_HC) - mean(D_SZ);
figure; cdfplot(abs(T)); title("Proposed Method"); xlabel("Absolute values of t-statistic");

proposed_common = p_bar;
for i = 1: length(proposed_common)
    proposed_common(i) = -log10(proposed_common(i));
end
for i = 1: 38 - length(proposed_common)
    proposed_common = [proposed_common, 0];
end
proposed_disc = p_tilde;
for i = 1: length(proposed_disc)
    proposed_disc(i) = -log10(proposed_disc(i));
end
for i = 1: 38 - length(proposed_disc)
    proposed_disc = [0, proposed_disc];
end

K = length(proposed_disc);

x = 1: K;
y1 = proposed_common;
y2 = proposed_disc;
bar_com = bar(x, y1);
hold on;
bar_disc = bar(x, y2);
p_line = yline(1.301,'--');

legend([bar_com, bar_disc, p_line],{'Common','Discriminative', 'p=0.05'}, 'Location', 'northwest');
xlabel("Component Index");
ylabel("-log_{10} p");

atoms_lambda_T_p = [atoms_lambda; T; [p_bar, p_tilde]];

load('../data/AOD1.mat');
cost_recon = norm([D_bar, D_tilde] * [Z_bar; Z_tilde] - Y', 'fro')^2/48546;
cost_D = cost_k(D, Y_range, K_bar, K_tilde, atoms_lambda);
cost_Z = best_hyperparams_list(1) * sum(sum(abs([Z_bar; Z_tilde])));

par2 = [150, 121];
class_labels_for_d = [ones(1, par2(1)), ones(1, par2(2)) + 1]; 
idx_a = 20;
idx_b = 32;
a = D(:, idx_a); b = D(:, idx_b);
[~, p_a, ~, ~] = ttest2(a(1:Y_range(2)), a(1 + Y_range(2):Y_range(3)),.05);

cost_a_if_common = 2 /(0.1 + exp(4 - bar_alpha * abs(atoms_lambda_T_p(2, idx_a)))) * h_k_cost_atom(a, class_labels_for_d, par2, 2);
cost_a_if_dis = 20 * exp(-tilde_alpha * abs(atoms_lambda_T_p(2, idx_a))) * f_k_cost_atom(a, class_labels_for_d, par2, 2);

cost_b_if_common = 2 /(0.1 + exp(4 - bar_alpha * abs(atoms_lambda_T_p(2, idx_b)))) * h_k_cost_atom(b, class_labels_for_d, par2, 2);
cost_b_if_dis = 20 * exp(-tilde_alpha * abs(atoms_lambda_T_p(2, idx_b))) * f_k_cost_atom(b, class_labels_for_d, par2, 2);


clear; close all;

function cost = f_k_cost_atom(d_k_atom_col, class_labels_for_d, M_per_class_list, eta_f)
    % Cost for one atom d_k to be discriminative.
    % Args:
    %   d_k_atom_col: A column vector representing the atom.
    %   class_labels_for_d: A vector of class labels for each sample.
    %   M_per_class_list: A vector of the number of samples per class.
    %   eta_f: A regularization parameter (default: 2.0).
    %   class1_val: Value of the first class (default: 0).

    m_k_overall = mean(d_k_atom_col);
    unique_classes = unique(class_labels_for_d);

    Sw_k = 0.0; % Within-class scatter for atom k
    Sb_k = 0.0; % Between-class scatter for atom k
    
    for c_idx = 1:length(unique_classes)
        C_val = unique_classes(c_idx);
        class_indices = find(class_labels_for_d == C_val);
        
        Mc = M_per_class_list(c_idx);
        d_k_c = d_k_atom_col(class_indices);
        m_c_k = mean(d_k_c);
        
        Sw_k = Sw_k + sum((d_k_c - m_c_k).^2);
        Sb_k = Sb_k + Mc * (m_c_k - m_k_overall)^2;
    end
    
    cost = Sw_k - Sb_k + eta_f * sum(d_k_atom_col.^2);
    fprintf("Fisher: Sw: %f, Sb: %f, norm: %f \n", Sw_k, Sb_k, eta_f * sum(d_k_atom_col.^2));
end

function cost = h_k_cost_atom(d_k_atom_col, class_labels_for_d, M_per_class_list, eta_h)
    % Cost for one atom d_k to be common.
    % Args:
    %   d_k_atom_col: A column vector representing the atom.
    %   class_labels_for_d: A vector of class labels for each sample.
    %   M_per_class_list: A vector of the number of samples per class.
    %   eta_h: A regularization parameter (default: 2.0).
    %   class1_val: Value of the first class (default: 0).

    m_k_overall = mean(d_k_atom_col);
    unique_classes = unique(class_labels_for_d);
    Sw_k = 0.0;
    Sb_k = 0.0;

    for c_idx = 1:length(unique_classes)
        C_val = unique_classes(c_idx);
        class_indices = find(class_labels_for_d == C_val);

        Mc = M_per_class_list(c_idx);
        d_k_c = d_k_atom_col(class_indices);
        m_c_k = mean(d_k_c);

        Sw_k = Sw_k + sum((d_k_c - m_c_k).^2);
        Sb_k = Sb_k + Mc * (m_c_k - m_k_overall)^2;
    end
    
    cost = Sb_k - Sw_k + eta_h * sum(d_k_atom_col.^2);
    fprintf("h: Sw: %f, Sb: %f, norm: %f \n", Sw_k, Sb_k, eta_h * sum(d_k_atom_col.^2));

end

function [cost_assigned] = cost_k(D, Y_range, K_bar, K_tilde, atoms_lambda)
    nClasses = 2;
    cost_assigned = [];
    for idx = 1: K_bar + K_tilde
        y_k = D(:, idx)';
        % [~, p, ~, ~] = ttest2(y_k(1:Y_range(2))', y_k(1 + Y_range(2):Y_range(3))',.05);
    
        m_k = mean(y_k);
        bt = 0; wn = 0;
        for c = 1: nClasses
            y_k_c = get_block_col(y_k, c, Y_range);
            Nc = size(y_k_c, 2);
            m_ck = mean(y_k_c);
            bt = bt + Nc * (m_ck - m_k)^2;
            for n = 1: Nc
                wn = wn + (y_k_c(n) - m_ck)^2/Nc;
            end
        end
        
        if idx <= K_bar
            cost = bt - wn + 2 * normF2(y_k);
            cost_assigned = [cost_assigned, atoms_lambda(idx) * cost];
        else
            cost = wn - bt + normF2(y_k);
            cost_assigned = [cost_assigned, atoms_lambda(idx) * cost];
        end
    end
end
function Mi = get_block_col(M, C, col_range)
	%% ================== File info ==========================
	% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
	% Time created	: 1/27/2016 2:39:32 AM
	% Last modified	: 1/27/2016 2:39:36 AM
	% Description	: Get column blocks of a big matrix M, blocks indexed by C 
	% 	INPUT
	%		M        : the big matrix. M = [M1 , M2 , ... , MC]
	%		C        : block indexes 
	%		col_range: a vector store the last index of each block. col_range(1) = 0.
	%					i-th block is indexed by col_range(i)+1: col_range(i+1).
	% 	OUTPUT 
	%		Mi: output block matrix  
	%
	%% ================== end File info ==========================

	id_sel = [];
	for i = 1: numel(C)
		c = C(i);
		id_sel = [id_sel, col_range(c) + 1: col_range(c+1)];
	end 
	Mi = M(:, id_sel, :);
end