function []= Plot_single()
%{Proposed Method %}
filePath = "..\experiments\real-permutate\bar_20_tilde_18\0.0027 0.175 0.34\2";
fileNames = ["D_bar.mat", "D_tilde.mat", "Z_bar.mat", "Z_tilde.mat", "info.mat"];
for i = 1: size(fileNames, 2)
    fn = fullfile(filePath, fileNames(i)); 
    load(fn);
end
Z = [Z_bar; Z_tilde];
D = [D_bar, D_tilde];
Sf = Z;
z_proposed = Z(13, :);
%{
JSTSP
%}
filePath = "..\data\ResultsofDLandICA";
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

z_sDL = Z(38, :);

R = corrcoef(z_proposed, z_sDL);

id = 13; % the id of map/component to plot
[~, p, ~, tstats] = ttest2(D(1:Y_range(2), id), D(1 + Y_range(2):Y_range(3), id),.05);
T = tstats.tstat;

load('indices.mat');

addpath('spm12');
fileName=['test']; %filename for saving the figures
numOfCV=size(Sf,1); %total number of estimated components

%% Select components of interest and z-score threshold
% Sf=-sR123.signal{3,2}(:,1:48546); %{dataset,run}(voxels) %% 1:48546___48547:97092___97093:145638
% Sf=ODD;

if ~exist('monthresh') %set z-score threshold (higher threshold means fewer plotted voxels
    monthresh = 2;
end

fprintf('Generating Data Display ... \n');

%% Create 4-D fMRI data and read in anatomical map
dat3f = zeros(numOfCV,53*63*46); 
dat3f(:,index_) = Sf(:,(1:length(index_)));
dat3f = reshape(dat3f,numOfCV,53,63,46);%why use reshape here?

fprintf('- Reading in Anatomical Data\n');
anat = spm_read_vols(spm_vol('nsingle_subj_T1_2_2_3.hdr'));
z_choice = 4:5:46;

tmp(:,:,:)=dat3f(id,:,:,z_choice);%x,y,time,slice?
figure;
xlabel('fMRI','FontSize', 12);
colorbar;
make_composite(sign(T)*tmp/stdN(tmp),anat(:,:,z_choice),monthresh); %plot component
% make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
clear tmp

%% save files as PDF
% 
% print('-depsc', fileName);
% hgsave(fileName);
% save(fileName);
% eps2pdf([fileName, '.eps'], 'C:\Program Files\gs\gs9.14\bin\gswin64.exe');   %for windows, the path of gs should be different
