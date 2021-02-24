
%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 2; % Number of UL UEs
J = 2; % Number of DL UEs
radar.codelength = 8;
%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
radar.TX = Mr;
radar.RX = Nr;
radar.noisepower = 0.01;
SNR.Btr = 2*ones(Nr,1);
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
SNR.rtr = 2*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.BB = 1;
SNR.BS_DL = 2*ones(J,1);
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
% Clutter 
SNR.CNR = ones(Nr,1);
radar.ell_max = 5; % algorithm 5

% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);

%% simulation parameters
[fdcomm_para, radar_fixed, radar_comm_para] = tsp_parameters(SNR,radar,fdcomm);
eta_rtr = radar_fixed.channelgain;
M = Mc + Mr;
ell_max = radar.ell_max;
dmode = 'radar';
SNR_rtr = -10:5:15;
radar_co_save = cell(length(SNR_rtr),1);
fdcomm_co_save = cell(length(SNR_rtr),1);
cov_co_save = cell(length(SNR_rtr),1);
radar_uc_save = cell(length(SNR_rtr),1);
fdcomm_uc_save = cell(length(SNR_rtr),1);
cov_uc_save = cell(length(SNR_rtr),1);
radar_rc_save = cell(length(SNR_rtr),1);
fdcomm_rc_save = cell(length(SNR_rtr),1);
cov_rc_save =  cell(length(SNR_rtr),1);
parfor ii = 1:length(SNR_rtr)
    %% Update the channel gain
    radar_comm_para_ii = radar_comm_para;
    radar_para = radar_fixed;
    eta_rtr = zeros(Mr,Nr);
    snrrtr = SNR_rtr(ii)*ones(Mr,Nr); 
    Sigma = zeros(M,M,Nr);
    eta_Btr = radar_comm_para_ii.Btrchannelgains;
    for nr = 1:Nr
        for mr = 1 : Mr
            eta_rtr(mr,nr) = db2pow(snrrtr(mr,nr))*radar_fixed.noisepower;
            Sigma(:,:,nr) = blkdiag(eta_rtr(:,nr).*eye(Mr),eta_Btr(nr)*eye(Mc));
        end
    end
    radar_para.channelgain = eta_rtr;
    radar_para.Sigma = Sigma;
    %% initialization 
    [fdcomm_ini, radar_ini, cov_temp] = tsp_ini_normal(radar_para,fdcomm_para,radar_comm_para_ii);
    %% WMMSE alternating algorithm
    fdcomm_co_op = fdcomm_ini;
    fdcomm_co = fdcomm_ini; 
    radar_co = radar_ini;
    radar_co_op = radar_ini;
    I_total_co = zeros(ell_max,1);
    I_total_co_op = zeros(ell_max,1);
    cov_temp_co = cov_temp;
    I_max_co = sum(fdcomm_co_op.MI_UL(:))+sum(fdcomm_co_op.MI_DL(:))+ sum(radar_co_op.MI_radar);
    % uncoded 
    radar_uc =  radar_ini;
    radar_uc_op = radar_uc;
    fdcomm_uc = fdcomm_ini;
    fdcomm_uc_op = fdcomm_ini;
    cov_temp_uc = cov_temp;
    I_total_uc = zeros(ell_max,1);
    I_total_uc_op = zeros(ell_max,1);
    I_max_uc = sum(fdcomm_uc_op.MI_UL(:))+sum(fdcomm_uc_op.MI_DL(:))+ sum(radar_uc.MI_radar);
    % random coded
    radar_rc =  radar_ini;
    radar_rc.codematrix = randomcode(radar_rc.Pr,Mr,radar_rc.codelength);
    radar_rc_op = radar_rc;
    fdcomm_rc = fdcomm_ini;
    fdcomm_rc_op = fdcomm_ini;
    cov_temp_rc = cov_temp;
    I_total_rc = zeros(ell_max,1);
    I_total_rc_op = zeros(ell_max,1);
    I_max_rc = sum(fdcomm_rc_op.MI_UL(:))+sum(fdcomm_rc_op.MI_DL(:))+ sum(radar_rc.MI_radar);
    ell = 1;
    while ell <= ell_max
        % Co-design WMMSE
        [fdcomm_co,radar_co_ast] = tsp_WMMSE_algorithm(fdcomm_co, radar_co, radar_comm_para_ii, cov_temp_co);
        % Update covariance matrices
        % Calculate A^star
        radar_co = tsp_Nearest_PAR(radar_co_ast);
        cov_temp_co = tsp_covmat(fdcomm_co,radar_co,radar_comm_para_ii);
        % Update linear receivers
        fdcomm_co = tsp_Comm_MMSE(fdcomm_co, radar_co, cov_temp_co);
        radar_co = tsp_radar_MMSE(radar_co, cov_temp_co);
        I_total_co(ell) = sum(fdcomm_co.MI_UL(:))+sum(fdcomm_co.MI_DL(:))+ sum(radar_co.MI_radar);
        if I_total_co(ell) > I_max_co
            I_total_co_op(ell) = I_total_co(ell);
            I_max_co = I_total_co(ell);
            radar_co_op = radar_co;
            fdcomm_co_op = fdcomm_co;
        else
            I_total_co_op(ell) = I_max_co;
        end
        fdcomm_co_op.I_total = I_total_co;
        fdcomm_co_op.I_total_op = I_total_co_op;
        cov_co_op = tsp_covmat(fdcomm_co_op,radar_co_op,radar_comm_para_ii);
        % uncoded WMMSE
        [fdcomm_uc, ~] =...
            tsp_WMMSE_algorithm_multimode(fdcomm_uc, radar_uc, radar_comm_para_ii, cov_temp_uc,dmode);
        cov_temp_uc = tsp_covmat(fdcomm_uc,radar_uc,radar_comm_para_ii);
        fdcomm_uc = tsp_Comm_MMSE(fdcomm_uc, radar_uc, cov_temp_uc);
        radar_uc = tsp_radar_MMSE(radar_uc, cov_temp_uc);
        I_total_uc(ell) = sum(fdcomm_uc.MI_UL(:))+sum(fdcomm_uc.MI_DL(:))+ sum(radar_uc.MI_radar);
        if I_total_uc(ell) > I_max_uc
            I_total_uc_op(ell) = I_total_uc(ell);
            I_max_uc = I_total_uc(ell);
            radar_uc_op = radar_uc;
            fdcomm_uc_op = fdcomm_uc;
        else
            I_total_uc_op(ell) = I_max_uc;
        end
        fdcomm_uc_op.I_total = I_total_uc;
        fdcomm_uc_op.I_total_op = I_total_uc_op;
        cov_uc_op = tsp_covmat(fdcomm_uc_op,radar_uc_op,radar_comm_para_ii);
        % random wmmse
        [fdcomm_rc, ~] =...
            tsp_WMMSE_algorithm_multimode(fdcomm_rc, radar_rc, radar_comm_para, cov_temp_rc,dmode);
        cov_temp_rc = tsp_covmat(fdcomm_rc,radar_rc,radar_comm_para);
        fdcomm_rc = tsp_Comm_MMSE(fdcomm_rc, radar_rc, cov_temp_rc);
        radar_rc = tsp_radar_MMSE(radar_rc, cov_temp_rc);
        I_total_rc(ell) = sum(fdcomm_rc.MI_UL(:))+sum(fdcomm_rc.MI_DL(:))+ sum(radar_rc.MI_radar);
        if I_total_rc(ell) > I_max_rc
            I_total_rc_op(ell) = I_total_rc(ell);
            I_max_rc = I_total_rc(ell);
            radar_rc_op= radar_rc;
            fdcomm_rc_op = fdcomm_rc;
        else
            I_total_rc_op(ell) = I_max_rc;
        end
        ell = ell+1;
        fdcomm_rc_op.I_total = I_total_rc;
        fdcomm_rc_op.I_total_op = I_total_rc_op;
        cov_rc_op = tsp_covmat(fdcomm_rc_op,radar_rc_op,radar_comm_para);
    end
    fdcomm_co_save{ii} = fdcomm_co_op;
    radar_co_save{ii} = radar_co_op;
    cov_co_save{ii} = cov_co_op;
    fdcomm_uc_save{ii} = fdcomm_uc_op;
    radar_uc_save{ii} = radar_uc_op;
    cov_uc_save{ii} = cov_uc_op;
    fdcomm_rc_save{ii} = fdcomm_rc_op;
    radar_rc_save{ii} = radar_rc_op;
    cov_rc_save{ii} = cov_rc_op;
    %file_name = ['radar_comparison', num2str(SNR_rtr(ii)),'dB.mat'];
   % save(file_name);
%    %savetofile('fdcomm_op_co','radar_op_co','cov_op_co',...
%        'fdcomm_op_uc','radar_op_uc','cov_op_uc','fdcomm_op_rc',...
%        'radar_op_rc','cov_op_rc',file_name)
%     
end
save('tsp_radar_snr_comparison.mat');
% function savetofile(data,fullfilename)
%     save(fullfilename,'data');
% end
 
