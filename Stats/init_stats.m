function [stats, RC_output] = init_stats()
stats=[];
stats.snr_asoc=[];%zeros(Nu,num_drops);
stats.LOS=[];%zeros(Nu*num_drops,Nb)
stats.snr_asoc_g2x=[];
stats.snr_asoc_x2u=[];
stats.pls=[];%zeros(Nu,num_drops);
stats.prx_asoc=[];%zeros(Nu,num_drops);
stats.cap=[];%zeros(Nu,num_drops);
stats.sinr=[];%zeros(Nu,num_drops);
stats.sinr_all=[];%zeros(Nu,num_iter,num_drops);
stats.snr=[];%zeros(Nu,Nb,num_drops);
stats.tput=[];
stats.idx_AS_ue=cell(0);%a cell array of size num_drops. Each entry is a cell array of size num_iter
stats.idx_AS_gnb=cell(0);
stats.cap_iter_drop=cell(0);%a cell array of size num_drops that holds matrices of size zeros(Nu,num_iter)
stats.macro2sc_inr=[];
stats.sc2macro_inr=[];
stats.nonSched_ues=[];

stats.scUEs_cap=[];
stats.scUEs_snr_asoc=[];
stats.scUEs_sinr=[];
stats.scUEs_pls=[];

stats.macroUEs_cap=[];
stats.macroUEs_snr_asoc=[];
stats.macroUEs_sinr=[];
stats.macroUEs_pls=[];

stats.rate = [];
stats.rate_rec = [];
stats.sinrdB_ave = [];
stats.sinrdB_periter = cell(0);
stats.mcs = cell(0);%a cell array of size num_drops. Each entry is a cell array of size num_iter
stats.rank = cell(0);%a cell array of size num_drops. Each entry is a cell array of size num_iter
stats.detect = cell(0);%a cell array of size num_drops. Each entry is a cell array of size num_iter

stats.UPT= []; % array to store the user perceived throughput of each UE
stats.latency = []; % array to store the average latency of each UE
stats.lambda = []; % average arrival rates for all UEs
stats.S = [] ;% average file sizes for all UEs
stats.success_rate = []; % array to store the successful transmission rate of each UE during each iteration

RC_output = [];
% if systemCfg.gNB.rate_control.DeltaUp>0
%     RC_output.Rate_rec = zeros(Nu,1);
%     RC_output.detect_outcome = nan(Nu, 1);
%     RC_output.Rate = zeros(Nu,1);
%     RC_output.SINR_lin_eff = zeros(Nu, 1);
%     RC_output.mcs_ues = nan(Nu, 1);
%     RC_output.rank_ues = zeros(1, Nu);
% end
