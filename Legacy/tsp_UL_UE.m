% Pd vs \vu for different SNR_Btr values

%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
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
fdcomm.DL_num = J;
SNR.BS_DL = 5*ones(J,1);
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.DL_num = J;
fdcomm.DLpower = J;
SNR.rtr = 5*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.BB = 1;
SNR.r_B = 2*ones(Mr,1);
radar.Pr = 1*ones(Mr,1);
radar.ell_max = 20; % algorithm 5
% Clutter 
SNR.CNR = ones(Nr,1);
SNR.r_DL = 1*ones(Mr,J);
II = 2:10; % # of UEs
M = 10;
radar_save = cell(length(II),M);
fdcomm_save = cell(length(II),M);
cov_save = cell(length(II),M);
for m = 1:M
    parfor ii = 1:length(II)
        I = II(ii);
        SNR_ii = SNR;
        fdcomm_ii = fdcomm;
        radar_ii = radar;
        SNR_ii.UL_DL = 1*ones(I,J);
        fdcomm_ii.UL_num = I;
        SNR_ii.UL_BS = 2*ones(I,1);
        SNR_ii.UL_r = 2*ones(I,Nr);
        SNR_ii.UL_BS = 2*ones(I,1);
        fdcomm_ii.ULpower = ones(I,1);
        %% priority weight equal
        fdcomm_ii.alpha_UL = 1/(Nr+I+J)*ones(I,1);
        fdcomm_ii.alpha_DL = 1/(Nr+I+J)*ones(J,1);
        radar_ii.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
        %% simulation parameters
        [fdcomm_para, radar_para, radar_comm_para] = tsp_parameters(SNR_ii,radar_ii,fdcomm_ii);
        [fdcomm_ini_ew_ni, radar_ini_ew_ni, cov_ew_ni] = tsp_ini_normal(radar_para,fdcomm_para,radar_comm_para);
        [fdcomm_ew_ni_op, radar_ew_ni_op, cov_ew_ni_op,~] =...
            tsp_altermating_projection_v1(fdcomm_ini_ew_ni, radar_ini_ew_ni,radar_comm_para,cov_ew_ni);
    %     I_total_op_ew_ni = fdcomm_ew_ni_op.I_total_op;
    %     I_UL_op_ew_ni = fdcomm_ew_ni_op.I_UL_op;
    %     I_DL_op_ew_ni = fdcomm_ew_ni_op.I_DL_op;
    %     I_radar_op_ew_ni = fdcomm_ew_ni_op.I_radar_op;
        radar_save{ii,m} = radar_ew_ni_op;
        fdcomm_save{ii,m} = fdcomm_ew_ni_op;
        cov_save{ii,m} = cov_ew_ni_op;
    end
end



 
