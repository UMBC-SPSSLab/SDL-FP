filePath = "..\data\ResultsofDLandICA";
% filePath = "../experiments\real\modified\16 24 65.9\runs\0.0033 0.1 0.35\1";
fileNames = ["ICAbestrun_Comp.mat", "ICAbestrun_act.mat"];
for i = 1: size(fileNames, 2)
    fn = fullfile(filePath, fileNames(i)); 
    load(fn);
end
D = A_hat_alldata;
Z = S_hat;
Y_range = [0, 121, 242];
K = size(D, 2);
M = size(D, 1);
V = size(Z, 2);

% Y_range = [0, 150, 271];

Sf = Z;
load('indices.mat');
addpath('spm12');%spm:statistical parametric mapping
fileName=['test']; %filename for saving the figures
monthresh = 2; %set z-score threshold (higher threshold means fewer plotted voxels
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
imagesc(D);
title("Imagesc of D");

% get the T values from two-sample t-test
pp = zeros(1, K);
idx = 1;
for k = 1: K
    [h(idx), pp(idx), ~, tstats(k)] = ttest2(D(1:Y_range(2),k), D(1 + Y_range(2):Y_range(3),k),.05);
    T(idx)=tstats(k).tstat;
    sele(idx) = sign(T(idx))*idx;
    idx = idx + 1;
end

% plot Z_bar
pp_mul = -100 + zeros(ceil(K/10), 10);
figure;
for s = 1: ceil(K/10)
    for j = 1: 10
        i = (s - 1) * 10 + j;
        if i <= K
            pp_mul(s, j) = pp(i);
            tmp(:, :, :) = dat3f(abs(sele(i)),:,:,z_choice);%x,y,time,slice?
            subplot_tight(ceil(K/10), 10, i, [0.001,0.005]);
            xlabel('fMRI', 'FontSize', 12);
            colorbar
            make_composite(sign(sele(i))*tmp/stdN(tmp),anat(:,:,z_choice),monthresh); %plot component
            % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
            clear tmp
        end
    end
end 
 
figure;
cdfplot(abs(Z(:)));
title("cdf of Z");

geomean_pp = geomean(pp, 'all');
geomedian_pp = geometric_median(pp);
fprintf("Geometric mean of p values of D: %f\n", geomean_pp);
fprintf("Geometric median of p values of D: %f\n", geomedian_pp);
disp("pp_mul");
disp(pp_mul);
clear; close all;