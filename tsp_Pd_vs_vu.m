%---------------- Simulation 2 experiment 1
%---------------- probability of detection vs radar
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
SNR_Btr = -15:5:5; 
radar.TX = Mr;
radar.RX = Nr;
radar.noisepower = 0.01;
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
SNR.rtr = 2*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
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
M = 5e3;
vu = 0:0.5:100;
PD = zeros(length(SNR_Btr),length(vu));
PFA = zeros(length(SNR_Btr),length(vu));
for ii = 1:length((SNR_Btr))
    snr_Btr = SNR_Btr(ii);
    SNR.Btr = snr_Btr*ones(Nr,1);
    for k = 1:length(vu)
        [~,~,radar_comm] = tsp_parameters(SNR,radar,fdcomm);
        file_name = ['uniformweight_BTR_', num2str(snr_Btr),'dB_co.mat'];
        load(file_name,'fdcomm','radar');
        [Pd,Pfa] = tsp_radar_detection(fdcomm,radar,radar_comm,M,vu(k));
        PD(ii,k) = Pd;
        PFA(ii,k) = Pfa;
    end
end