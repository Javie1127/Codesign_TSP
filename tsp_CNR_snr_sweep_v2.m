% Pd vs \vu for different SNR_Btr values

%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas
radar.TX = Mr;
radar.RX = Nr;      
radar.noisepower = 0.01;
K = 8;
radar.codelength = K;
%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 2; % Number of UL UEs
J = 2; % Number of DL UEs
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
%% Set the SNRs in dB
%SNR.rtr = randi([-5,5],Mr,Nr);
SNR.rtr = 5*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.Btr = 1*ones(Nr,1); 
SNR.BB = 1;
SNR.BS_DL = 2*ones(J,1);
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
radar.ell_max = 10; % algorithm 5
% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
%% simulation parameters
[fdcomm_para, radar_para, radar_comm_para] = tsp_parameters_no_clutter(SNR,radar,fdcomm);
eta_rtr = radar_para.channelgain;
ell_max = radar.ell_max;
%% Parameters to save
CNR_sweep = -30:5:30;
radar_save = cell(length(CNR_sweep),1);
fdcomm_save = cell(length(CNR_sweep),1);
cov_save = cell(length(CNR_sweep),1);
%% Loop through CNRs 
Sigma_c = db2pow(CNR_sweep)*sigma_0;
parfor ii = 1:length(CNR_sweep)
    radar_para_ii = radar_para;
    fdcomm_para_ii = fdcomm_para;
    radar_comm_para_ii = radar_comm_para;
    %% Update the Clutter cov matrix
    %CNR = CNR_sweep(ii);
    R_C = zeros(K,K,Nr); 
    rho = 0.5;
    % Method 1 C is a gaussian rv (0,sigma_c);
    %sigma_c = db2pow(CNR)*sigma_0;
    sigma_c = Sigma_c(ii);
    for nr = 1:Nr
        R_C_nr = zeros(K,K);
        for kk = 1 : K 
            for jj = 1:K
                R_C_nr(kk,jj) = sigma_c*rho^(abs(kk-jj));
            end
        end
        R_C(:,:,nr) = R_C_nr;
    end
    radar_para_ii.cluttercov = R_C;
    % radar TX - BS 
    mu_r_BS = 0.5*ones(Mr,1);
    kappa_BS = 0.5;
    H_r_BS = zeros(Nc,Mr);
    eta_rB = zeros(Mr,1);
    for mr = 1 : Mr 
        eta_rB(mr) = 0.8*sigma_c;
        H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
    end
    radar_comm_para_ii.radar2BSchannels = H_r_BS;
    radar_comm_para_ii.r2Bchannelgains = eta_rB;
    mu_DL = 0.2*ones(Mr,J);
    kappa_DL = 0.3;
    H_r_DL = cell(J,1);
    eta_r_DL = zeros(Mr,J);
    N_DL = 2*ones(J,1); % Number of the DL UE antennas
    for jj = 1 : J
        H_r_j = zeros(N_DL(jj),Mr);
        for mr = 1:Mr
            eta_r_DL(mr,jj) = 0.7*sigma_c; 
            H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
        end
        H_r_DL{jj,1} = H_r_j;
    end
    radar_comm_para_ii.radar2DLchannnels = H_r_DL;
    radar_comm_para_ii.radar2DLchannelgains = eta_r_DL;
    %% initialization 
    [fdcomm_ini, radar_ini, cov_temp] = tsp_ini_normal(radar_para_ii,fdcomm_para_ii,radar_comm_para_ii);
    [fdcomm_op, radar_op,cov_op, ~] = tsp_altermating_projection_v1(fdcomm_ini, radar_ini,radar_comm_para_ii,cov_temp);
    radar_save{ii} = radar_op;
    fdcomm_save{ii} = fdcomm_op;
    cov_save{ii} = cov_op;
end

 
