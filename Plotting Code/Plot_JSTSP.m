filePath = "..\data\ResultsofDLandICA";
% filePath = "../experiments\real\modified\16 24 65.9\runs\0.0033 0.1 0.35\1";
fileNames = ["DLbestrun_Common_Act.mat", "DLbestrun_Common_Comp.mat", "DLbestrun_Disc_Act.mat", "DLbestrun_Disc_Comp.mat"];
for i = 1: size(fileNames, 2)
    fn = fullfile(filePath, fileNames(i)); 
    load(fn);
end
D_bar = D0; D_tilde = D; Z_bar = X0; Z_tilde = X;
Y_range = [0, 121, 242];
[K_bar, K_tilde] = deal(size(D_bar, 2), size(D_tilde, 2));
K = K_bar + K_tilde;
M = size(D_bar, 1);
V = size(D_bar, 2);

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
dat3f = reshape(dat3f, K, 53, 63, 46);%why use reshape here?

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
idx = 1;
for k = 1: K_bar
    [h(idx), p_bar(idx), ~, tstats(k)] = ttest2(D_bar(1:Y_range(2),k), D_bar(1 + Y_range(2):Y_range(3),k),.05);
    T(idx)=tstats(k).tstat;
    sele(idx) = sign(T(idx))*idx;
    idx = idx + 1;
end

for k = 1: K_tilde
    [h(idx), p_tilde(idx - K_bar), ~, tstats(k)] = ttest2(D_tilde(1:Y_range(2),k), D_tilde(1 + Y_range(2):Y_range(3),k),.05);
    T(idx) = tstats(k).tstat;
    sele(idx) = sign(T(idx))*idx;
    idx = idx + 1;
end

% plot Z_bar
pp_bar = -100 + zeros(ceil(K_bar/10), 10);
figure;
for s = 1: ceil(K_bar/10)
    for j = 1: 10
        i = (s - 1) * 10 + j;
        if i <= K_bar
            pp_bar(s, j) = p_bar(i);
            tmp(:, :, :) = dat3f(abs(sele(i)),:,:,z_choice);%x,y,time,slice?
            subplot_tight(ceil(K_bar/10), 10, i, [0.001,0.005]);
            xlabel('fMRI', 'FontSize', 12);
            colorbar
            make_composite(sign(sele(i))*tmp/stdN(tmp),anat(:,:,z_choice),monthresh); %plot component
            % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
            clear tmp
        end
    end
end 
 
% plot Z_tilde
pp_tilde = -100 + zeros(ceil(K_tilde/10), 10);
figure;
for s = 1: ceil(K_tilde/10)
    for j = 1: 10
        i = (s - 1) * 10 + j;
        if i <= K_tilde
            pp_tilde(s, j) = p_tilde(i);
            tmp(:,:,:) = dat3f(abs(sele(K_bar + i)), :, :, z_choice);%x,y,time,slice?
            subplot_tight(ceil(K_tilde/10), 10, i, [0.001,0.005]);
            xlabel('fMRI', 'FontSize', 12);
            colorbar
            make_composite(sign(sele(i + K_bar))*tmp/stdN(tmp), anat(:, :, z_choice), monthresh); %plot component
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
figure; cdfplot(abs(T)); title("Supervised DL"); xlabel("Absolute values of t-statistic");

 clear; close all;