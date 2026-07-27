
% clear
%load('./results/LRSDL_mod_val/AOD1SIRP23fusion_N_60_k_16_k0_1_l1_2.5e-06_l2_0.17_0126_132609.mat')

%% P1
% num_comp=size(X,1);
% for k=1:num_comp
% %     [h(k),p(k),~,tstats(k)] = ttest2(X(k,1:60),X(k,61:120),.05);
%      [h(k),p(k),~,tstats(k)] = ttest2(X(k,1:30),X(k,31:60),.05);
%     T(k)=tstats(k).tstat;
% end
% T
% p
% 
% num_comp=size(X0,1);
% for k=1:num_comp
% %     [h0(k),p0(k),~,tstats(k)] = ttest2(X0(k,1:60),X0(k,61:120),.05);
%     [h0(k),p0(k),~,tstats(k)] = ttest2(X0(k,1:30),X0(k,31:60),.05);
%     T0(k)=tstats(k).tstat;
% end
% T0
% p0

%% P2
num_comp=size(D,2);
for k=1:num_comp
    [h(k),p(k),~,tstats(k)] = ttest2(D(1:121,k),D(122:242,k),.05);
    T(k)=tstats(k).tstat;
end
T
p

% num_comp=size(D0,2);
% for k=1:num_comp
%     [h0(k),p0(k),~,tstats(k)] = ttest2(D0(1:121,k),D0(122:242,k),.05);
%     T0(k)=tstats(k).tstat;
% end
% T0
% p0

%% ICA
% 
% num_comp=size(best_ica_test_activations,2);
% % param
% for k=1:num_comp
%     [h(k),p(k),~,tstats(k)] = ttest2(best_ica_test_activations(1:45,k),best_ica_test_activations(46:81,k),.05);
%     T(k)=tstats(k).tstat;
% end
% disp('test t statistic')
% T
% p
% 
% %% ICA train
% num_comp=size(best_ica_train_activations,2);
% % param
% for k=1:num_comp
%     [h(k),p(k),~,tstats(k)] = ttest2(best_ica_train_activations(1:105,k),best_ica_train_activations(106:end,k),.05);
%     T(k)=tstats(k).tstat;
% end
% 
% disp('train t statistic')
% T
% p
% 
% 

