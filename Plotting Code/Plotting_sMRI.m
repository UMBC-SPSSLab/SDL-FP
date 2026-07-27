clear
clc

%% Load results and true voxel locations (indices)
load('C:\Users\Yuri Levin-Schwartz\Desktop\Fall 2016\IVA_tIVA\_Results123_IVA_new_Fin.mat') %load results
load('C:\Users\Yuri Levin-Schwartz\Desktop\Fall 2015\Order 24\inDEx.mat') % load indices
ind_s=find(index_true); %determine indices to plot

fileName=['test']; %filename for saving the figures
numOfCV=24; %total number of estimated components

%% Select components of interest and z-score threshold
Ss=-sR123.signal{3,2}(:,1:48546); %{dataset,run}(voxels) %% 1:48546___48547:97092___97093:145638

if ~exist('monthresh') %set z-score threshold (higher threshold means fewer plotted voxels
    monthresh = 2.7
end

fprintf('Generating Data Display ... \n');
if ~exist('sele') %set the components to be plotted
    sele = [5 6 12 15 19]
end

%% Create 4-D sMRI data and read in anatomical map
dat3s = zeros(numOfCV,91*109*91);
dat3s(:,ind_s) = Ss(:,(1:length(ind_s)));
dat3s = reshape(dat3s,numOfCV,91,109,91);

fprintf('- Reading in Anatomical Data\n');
anat = spm_read_vols(spm_vol('C:\Users\Yuri Levin-Schwartz\Desktop\Spring 2013\Signal Processing\Fusion Project\From the Server\nsingle_subj_T1_2_2_3.img'));

% select slices to show ("1:1:46" will show all 46 slices)
z_choice = 4:5:46;

%% Reslice sMRI so that it is in the same space as the fMRI

Vf = spm_vol('C:\Users\Yuri Levin-Schwartz\Desktop\con_0009.img');
% Load structual MRI header to create temp files
% - structural maps based on projections are same size and shape as
%   original structural data
V  = spm_vol('C:\Users\Yuri Levin-Schwartz\Desktop\swc1T1_0001.nii');

% Cycle through the components
for ind = 1:numOfCV
    outfilename = ['C:\Users\Yuri Levin-Schwartz\Desktop\structtemp\sMRI_tmp' num2str(ind) '.img'];
    icatb_write_nifti_data(outfilename,V,squeeze(dat3s(ind,:,:,:)));
    VS = spm_vol(outfilename);    
    spm_reslice([Vf VS]);
    
    % Read in reslice data
    Vr = spm_vol(['C:\Users\Yuri Levin-Schwartz\Desktop\structtemp\rsMRI_tmp' num2str(ind) '.img']);
    stemp = spm_read_vols(Vr);
    ndat3s(ind,:,:,:) = stemp;
    
end

%% Plot individual images
for i=1:length(sele) %1:numOfCV
    
    tmp2(:,:,:)=ndat3s(abs(sele(i)),:,:,z_choice);  

    figure
    make_composite(sign(sele(1))*tmp2/stdN(tmp2),anat(:,:,z_choice),monthresh);
    clear tmp2
end
%% save files as PDF
% 
% print('-depsc', fileName);
% hgsave(fileName);
% save(fileName);
% eps2pdf([fileName, '.eps'], 'C:\Program Files\gs\gs9.14\bin\gswin64.exe');   %for windows, the path of gs should be different
