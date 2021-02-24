%%% Convergence analysis

%% Array Parameters
Mr = 6; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 4; % Number of UL UEs
J = 2; % Number of DL UEs

%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
SNR_rtr = -10:5:15;
radar.TX = Mr;
radar.RX = Nr;
radar.noisepower = 0.01;
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
SNR.Bmr = 1*ones(Nr,1);
SNR.Btr = 1*ones(Nr,1); 
SNR.BB = 1;
SNR.BS_DL = 5*ones(J,1);
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
% Clutter 
SNR.CNR = ones(Nr,1);
radar.ell_max = 10; % algorithm 5

% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
for snr_rtr = SNR_rtr
    SNR.rtr = snr_rtr*ones(Mr,Nr);
    %% WMMSE alternating algorithm
    [fdcomm, radar] = tsp_altermating_projection(SNR, fdcomm, radar,'normal_ini');
    file_name = ['uniformweight_radar_', num2str(snr_rtr),'dB_.mat'];
    save(file_name);
end

 
