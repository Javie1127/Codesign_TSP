%---------------- Simulation 2 experiment 2
%---------------- probability of detection vs radar snr
load('tsp_radar_snr_comparison_v5.mat')
M = 1e4;
vu = 0:0.5:100;
PD_co = zeros(length(SNR_rtr),length(vu));
PFA_co = zeros(length(SNR_rtr),length(vu));
PD_co_wc = zeros(length(SNR_rtr),length(vu));
PFA_co_wc = zeros(length(SNR_rtr),length(vu));
PD_uc = zeros(length(SNR_rtr),length(vu));
PFA_uc = zeros(length(SNR_rtr),length(vu));
PD_uc_wc = zeros(length(SNR_rtr),length(vu));
PFA_uc_wc = zeros(length(SNR_rtr),length(vu));
PD_rc = zeros(length(SNR_rtr),length(vu));
PFA_rc = zeros(length(SNR_rtr),length(vu));
PD_rc_wc = zeros(length(SNR_rtr),length(vu));
PFA_rc_wc = zeros(length(SNR_rtr),length(vu));
radar_comm_para_co_wc = radar_comm_para_co;
radar_comm_para_co_wc.Btrchannelgains = zeros(Nr,1);
for ii = 1:length((SNR_rtr))
    fdcomm_co_ii = fdcomm_co_save{ii};
    radar_co_ii = radar_co_save{ii};
    cov_co_ii = cov_co_save{ii};
    fdcomm_co_wc_ii = fdcomm_co_wc_save{ii};
    radar_co_wc_ii = radar_co_wc_save{ii};
    cov_co_wc_ii = cov_co_wc_save{ii};
    fdcomm_uc_ii = fdcomm_uc_save{ii};
    radar_uc_ii = radar_uc_save{ii};
    cov_uc_ii = cov_uc_save{ii};
    fdcomm_uc_wc_ii = fdcomm_uc_wc_save{ii};
    cov_uc_wc_ii = cov_uc_wc_save{ii};
    fdcomm_rc_ii = fdcomm_rc_save{ii};
    radar_rc_ii = radar_rc_save{ii};
    cov_rc_ii = cov_rc_save{ii};
    fdcomm_rc_wc_ii = fdcomm_rc_wc_save{ii};
    radar_rc__wc_ii = radar_rc_wc_save{ii};
    cov_rc_wc_ii = cov_rc_wc_save{ii};
    PD_co_ii = zeros(length(vu),1);
    PFA_co_ii = zeros(length(vu),1);
    PD_co_wc_ii = zeros(length(vu),1);
    PFA_co_wc_ii = zeros(length(vu),1);
    PD_uc_ii = zeros(length(vu),1);
    PFA_uc_ii = zeros(length(vu),1);
    PD_uc_wc_ii = zeros(length(vu),1);
    PFA_uc_wc_ii = zeros(length(vu),1);
    PD_rc_ii = zeros(length(vu),1);
    PFA_rc_ii = zeros(length(vu),1);
    PD_rc_wc_ii = zeros(length(vu),1);
    PFA_rc_wc_ii = zeros(length(vu),1);
   parfor k = 1:length(vu)
        [Pd_co,Pfa_co] = tsp_radar_detection(fdcomm_co_ii,radar_co_ii,...
            radar_comm_para_co,M,vu(k),cov_co_ii);
        PD_co_ii(k) = Pd_co;
        PFA_co_ii(k) = Pfa_co;
        [Pd_co_wc,Pfa_co_wc] = tsp_radar_detection(fdcomm_co_ii,radar_co_ii,...
            radar_comm_para_co_wc,M,vu(k),cov_co_ii);
        PD_co_wc_ii(k) = Pd_co_wc;
        PFA_co_wc_ii(k) = Pfa_co_wc;
        %% uniformly coded 
        [Pd_uc,Pfa_uc] = tsp_radar_detection(fdcomm_uc_ii,radar_uc_ii,...
            radar_comm_para_co,M,vu(k),cov_uc_ii);
        PD_uc_ii(k) = Pd_uc;
        PFA_uc_ii(k) = Pfa_uc;
        [Pd_uc_wc,Pfa_uc_wc] = tsp_radar_detection(fdcomm_uc_ii,radar_uc_ii,...
            radar_comm_para_co_wc,M,vu(k),cov_uc_ii);
        PD_uc_wc_ii(k) = Pd_uc_wc;
        PFA_uc_wc_ii(k) = Pfa_uc_wc;
        %% randomly coded 
        [Pd_rc, Pfa_rc] = tsp_radar_detection(fdcomm_rc_ii,radar_rc_ii,...
            radar_comm_para_co,M,vu(k),cov_rc_ii);
        PD_rc_ii(k) = Pd_rc;
        PFA_rc_ii(k) = Pfa_rc;
        [Pd_rc_wc, Pfa_rc_wc] = tsp_radar_detection(fdcomm_rc_ii,radar_rc_ii,...
            radar_comm_para_co_wc,M,vu(k),cov_rc_ii);
        PD_rc_wc_ii(k) = Pd_rc_wc;
        PFA_rc_wc_ii(k) = Pfa_rc_wc;
    end
    PD_co(ii,:) = PD_co_ii.';
    PD_co_wc(ii,:) = PD_co_wc_ii.';
    PFA_co(ii,:) = PFA_co_ii.';
    PFA_co_wc(ii,:) = PFA_co_wc_ii.';
    PD_uc(ii,:) = PD_uc_ii.';
    PD_uc_wc(ii,:) = PD_uc_wc_ii;
    PFA_uc(ii,:) = PFA_uc_ii.';
    PFA_uc_wc(ii,:) = PFA_uc_wc_ii.';
    PD_rc(ii,:) = PD_rc_ii.';
    PD_rc_wc(ii,:) = PD_rc_wc_ii.';
    PFA_rc(ii,:) = PFA_rc_ii.';
    PFA_rc_wc(ii,:) = PFA_rc_wc_ii.';
end
save('tsp_Pd_vs_radar_SNR_v3.mat')