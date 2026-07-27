clear 
clc

%% load ERP components and original EEG data and specify the components to display

load('/media/krishna/OS/acads/5_fmri_project/lrsdl/data/Task_data/AOD1.mat'); % load ERP results
% load('C:\Users\Yuri Levin-Schwartz\Desktop\Corrected fMRI and EEG\EEG_data_corrected.mat'); % load orignal EEG data

Se=ODD;
if ~exist('sele') %set the components to be plotted
    sele = [6 7]
end

%% Plot ERP component (using timepoints 50:500)
eind=-114:2:786;

for i=1:length(sele) %1:numOfCV
    
    tmp = detrend(Se(abs(sele(i)),:)); %Detrend the ERP components
    figure
    hp = plot(eind,sign(sele(i))*tmp/stdN(tmp),'.-r');
    set(hp,'linewidth',2);
    hl = legend('ERP Component', 'Location', 'NorthOutside');
    set(hl,'Interpreter','none')
    xlabel('milliseconds','FontSize', 12);
    ylabel('microvolts','FontSize', 12);

    axis ij; 
    axis([-115.    788.   -4    4]);
    set(gca,'FontSize', 12);
    clear tmp;
end

%% save files as PDF
% 
% print('-depsc', fileName);
% hgsave(fileName);
% save(fileName);
% eps2pdf([fileName, '.eps'], 'C:\Program Files\gs\gs9.14\bin\gswin64.exe');   %for windows, the path of gs should be different