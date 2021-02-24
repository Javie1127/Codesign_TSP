%---------------- Simulation 2 experiment 2
%---------------- probability of detection vs radar snr
load('tsp_radar_snr_comparison_v3.mat')
M = 4e3;
vu = 0:0.2:500;
PD_co = zeros(length(SNR_rtr),length(vu));
PFA_co = zeros(length(SNR_rtr),length(vu));
PD_uc = zeros(length(SNR_rtr),length(vu));
PFA_uc = zeros(length(SNR_rtr),length(vu));
PD_rc = zeros(length(SNR_rtr),length(vu));
PFA_rc = zeros(length(SNR_rtr),length(vu));
for ii = 1:length((SNR_rtr))
    fdcomm_co_ii = fdcomm_co_save{ii};
    radar_co_ii = radar_co_save{ii};
    cov_co_ii = cov_co_save{ii};
    fdcomm_uc_ii = fdcomm_uc_save{ii};
    radar_uc_ii = radar_uc_save{ii};
    cov_uc_ii = cov_uc_save{ii};
    fdcomm_rc_ii = fdcomm_rc_save{ii};
    radar_rc_ii = radar_rc_save{ii};
    cov_rc_ii = cov_rc_save{ii};
    PD_co_ii = zeros(length(vu),1);
    PFA_co_ii = zeros(length(vu),1);
    PD_uc_ii = zeros(length(vu),1);
    PFA_uc_ii = zeros(length(vu),1);
    PD_rc_ii = zeros(length(vu),1);
    PFA_rc_ii = zeros(length(vu),1);
    for k = 1:length(vu)
        [Pd_co,Pfa_co] = tsp_radar_detection(fdcomm_co_ii,radar_co_ii,...
            radar_comm_para,M,vu(k),cov_co_ii);
        PD_co_ii(k) = Pd_co;
        PFA_co_ii(k) = Pfa_co;
        [Pd_uc,Pfa_uc] = tsp_radar_detection(fdcomm_uc_ii,radar_uc_ii,...
            radar_comm_para,M,vu(k),cov_uc_ii);
        PD_uc_ii(k) = Pd_uc;
        PFA_uc_ii(k) = Pfa_uc;
        [Pd_rc, Pfa_rc] = tsp_radar_detection(fdcomm_rc_ii,radar_rc_ii,...
            radar_comm_para,M,vu(k),cov_rc_ii);
        PD_rc_ii(k) = Pd_rc;
        PFA_rc_ii(k) = Pfa_rc;
    end
    PD_co(ii,:) = PD_co_ii.';
    PFA_co(ii,:) = PFA_co_ii.';
    PD_uc(ii,:) = PD_uc_ii.';
    PFA_uc(ii,:) = PFA_uc_ii.';
    PD_rc(ii,:) = PD_rc_ii.';
    PFA_rc(ii,:) = PFA_rc_ii.';
end
save('tsp_Pd_vs_radar_SNR_v1.mat')