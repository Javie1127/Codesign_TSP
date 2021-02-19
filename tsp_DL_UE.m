% Pd vs \vu for different SNR_Btr values

%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 4; % Number of UL UEs
J = 2; % Number of DL UEs
K = 8;
radar.codelength = K;
%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
SNR.Btr = 2*ones(Nr,1);
radar.TX = Mr;
radar.RX = Nr;
sigma_0 = 0.01;
radar.noisepower = sigma_0;
SNR.Btr = ones(Nr,1);
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
SNR.rtr = 2*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.BB = 1;
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
radar.Pr = 1*ones(Mr,1);
radar.ell_max = 20; % algorithm 5
% Clutter 
SNR.CNR = ones(Nr,1);
JJ = 2:10; % # of UEs
M = 10;
radar_save = cell(length(JJ),M);
fdcomm_save = cell(length(JJ),M);
cov_save = cell(length(JJ),M);
for m = 1:M
    parfor jj = 1:length(JJ)
        J = JJ(jj);
        fdcomm_jj = fdcomm;
        SNR_jj = SNR;
        radar_jj = radar;
        fdcomm_jj.DL_num = J;
        SNR_jj.BS_DL = 5*ones(J,1);
        SNR_jj.UL_DL = 1*ones(I,J);
        SNR_jj.r_DL = 1*ones(Mr,J);
        fdcomm_jj.DLpower = J;
        %% priority weight equal
        fdcomm_jj.alpha_UL = 1/(Nr+I+J)*ones(I,1);
        fdcomm_jj.alpha_DL = 1/(Nr+I+J)*ones(J,1);
        radar_jj.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
        %% simulation parameters
        [fdcomm_para, radar_para, radar_comm_para] = tsp_parameters(SNR_jj,radar_jj,fdcomm_jj);
        [fdcomm_ini_ew_ni, radar_ini_ew_ni, cov_ew_ni] = tsp_ini_normal(radar_para,fdcomm_para,radar_comm_para);
        [fdcomm_ew_ni_op, radar_ew_ni_op, cov_ew_ni_op,~] =...
            tsp_altermating_projection_v1(fdcomm_ini_ew_ni, radar_ini_ew_ni,radar_comm_para,cov_ew_ni);
    %     I_total_op_ew_ni = fdcomm_ew_ni_op.I_total_op;
    %     I_UL_op_ew_ni = fdcomm_ew_ni_op.I_UL_op;
    %     I_DL_op_ew_ni = fdcomm_ew_ni_op.I_DL_op;
    %     I_radar_op_ew_ni = fdcomm_ew_ni_op.I_radar_op;
        radar_save{jj,m} = radar_ew_ni_op;
        fdcomm_save{jj,m} = fdcomm_ew_ni_op;
        cov_save{jj,m} = cov_ew_ni_op;
    end
end



 
