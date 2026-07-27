function []= Plot_fMRI(X,sele,monthresh)
% clear
% clc
% close all
%% Load results and true voxel locations (indices)
% load('/media/krishna/OS/acads/5_fmri_project/lrsdl/data/Task_data/AOD3.mat') %load results
% load('D:\Research and Study\Materials and Codes\RP\fMRI\data\asilomar\code\Task_data\indices.mat') %load results
% load('/media/krishna/OS/acads/5_fmri_project/lrsdl/results/AOD1/k30ko20lam1_0.01lam2_0.01lam3_0.mat')
% load('/media/krishna/OS/acads/5_fmri_project/lrsdl/data/Task_data/fMRI_data_corrected.mat') % load indices
Sf=X;
load('indices.mat');
%load('indices.mat');
addpath('spm12');
fileName=['test']; %filename for saving the figures
numOfCV=size(Sf,1); %total number of estimated components

%% Select components of interest and z-score threshold
% Sf=-sR123.signal{3,2}(:,1:48546); %{dataset,run}(voxels) %% 1:48546___48547:97092___97093:145638
% Sf=ODD;

if ~exist('monthresh') %set z-score threshold (higher threshold means fewer plotted voxels
    monthresh = 2
end

fprintf('Generating Data Display ... \n');
if ~exist('sele') %set the components to be plotted
    sele = [1 3 4]
end

%% Create 4-D fMRI data and read in anatomical map
dat3f = zeros(numOfCV,53*63*46); 
dat3f(:,index_) = Sf(:,(1:length(index_)));
dat3f = reshape(dat3f,numOfCV,53,63,46);%why use reshape here?
% dat3f = reshape(dat3f,numOfCV,54,29,31);

fprintf('- Reading in Anatomical Data\n');
anat = spm_read_vols(spm_vol('nsingle_subj_T1_2_2_3.hdr'));
%anat is the graph or image of anatomical sturcture?
% select slices to show ("1:1:46" will show all 46 slices)
z_choice = 4:5:46;

%for individual images
for i=1:length(sele) %1:numOfCV
    tmp(:,:,:)=dat3f(abs(sele(i)),:,:,z_choice);%x,y,time,slice?
%     m=(sign(sele(i))*Sf/stdN(Sf)); % use if you want to normalize based on all voxels (as opposed to only the plotted ones (see line 43))

    figure;
    xlabel('fMRI','FontSize', 12);
    colorbar
    make_composite(sign(sele(i))*tmp/stdN(tmp),anat(:,:,z_choice),monthresh); %plot component
    % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
    clear tmp
end

%% save files as PDF
% 
% print('-depsc', fileName);
% hgsave(fileName);
% save(fileName);
% eps2pdf([fileName, '.eps'], 'C:\Program Files\gs\gs9.14\bin\gswin64.exe');   %for windows, the path of gs should be different
