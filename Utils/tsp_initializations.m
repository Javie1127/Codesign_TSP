function [fdcomm, radar,radar_comm,cov] = tsp_initializations(systemCfg)
%% Radar-only parameters
radar.Tx =systemCfg.Antennas.Num_Radar_Txs;
radar.Rx =systemCfg.Antennas.Num_Radar_Rxs;
radar.codelength = systemCfg.System_Parameters.Num_PRI;
radar.num_range_cell = systemCfg.System_Parameters.Num_Range_Cell;
radar.gamma_r = systemCfg.Powers.radar_PAR; % PAR level
radar.K_factor = systemCfg.Channel_Modeling.K_Rician_radar;
radar.Rician_direct = systemCfg.Channel_Modeling.Rician_direct_radar;
radar.SNR = systemCfg.Powers.SNR_radar;
radar.CNR = systemCfg.Powers.CNR;
radar.PAR = systemCfg.Powers.radar_PAR;
%% FD-only parameters
fdcomm.BSTx = systemCfg.Antennas.Num_BS_Antennas;
fdcomm.BSRx = systemCfg.Antennas.Num_BS_Antennas;
fdcomm.UL_SNR = systemCfg.Powers.SNR_UL;
fdcomm.DL_SNR = systemCfg.Powers.SNR_DL;
fdcomm.Weights = systemCfg.Cases.Weights;
fdcomm.Initializations = systemCfg.Cases.Initializations;
if systemCfg.Cases.FD.isEnabled
    fdcomm.UL_num = systemCfg.Cases.FD.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.FD.Num_DL_UE;
    fdcomm.UL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas*ones(fdcomm.UL_num,1);
    fdcomm.DL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas*ones(fdcomm.DL_num,1);
    fdcomm.ULstream_num = fdcomm.UL_UE_Ant;
    fdcomm.DLstream_num = fdcomm.DL_UE_Ant;
    radar_comm.isCollaborate = systemCfg.Cases.FD.isCollaborationEnabled;
    fdcomm.theta_BT = systemCfg.Channel_Modeling.theta_Bt;
elseif systemCfg.Cases.DL.isEnabled
    fdcomm.UL_num = systemCfg.Cases.DL.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.DL.Num_DL_UE;
    fdcomm.DL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas*ones(fdcomm.DL_num,1);
    fdcomm.UL_UE_Ant = 0;
    fdcomm.DLstream_num = fdcomm.DL_UE_Ant;
    radar_comm.isCollaborate = systemCfg.Cases.DL.isCollaborationEnabled;
    fdcomm.theta_BT = systemCfg.Channel_Modeling.theta_Bt;
elseif systemCfg.Cases.UL.isEnabled
    fdcomm.UL_num = systemCfg.Cases.UL.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.UL.Num_DL_UE;
    fdcomm.DL_UE_Ant = 0;
    fdcomm.UL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas*ones(fdcomm.UL_num,1);
    fdcomm.ULstream_num = fdcomm.UL_UE_Ant;
    radar_comm.isCollaborate = 0;
end
fdcomm.pathloss = systemCfg.Channel_Modeling.Path_Loss;
fdcomm.K_factor = systemCfg.Channel_Modeling.K_Rician_SI;
fdcomm.UL_pathloss = systemCfg.Powers.SNR_UL;
% noise power
Lp = systemCfg.System_Parameters.Noise_density;
Bw = systemCfg.System_Parameters.System_Bandwidth;
NF_BS = systemCfg.System_Parameters.BS_Noise_Figure;
F_BS = 10^((NF_BS)/10);
fdcomm.BS_noise_power = F_BS*0.001*10^(Lp/10)*Bw;
fdcomm.BS_power = 10^(fdcomm.DL_SNR/10)*fdcomm.BS_noise_power;

NF_UE = systemCfg.System_Parameters.UE_Noise_Figure;
F_UE = 10^((NF_UE)/10);
fdcomm.UE_noise_power = F_UE*0.001*10^(Lp/10)*Bw;
fdcomm.UL_power = 10^(fdcomm.UL_SNR/10)*fdcomm.BS_noise_power*ones(fdcomm.UL_num,1); % per UE power

NF_radar = systemCfg.System_Parameters.Radar_Noise_Figure;
F_radar = 10^((NF_radar)/10);
radar.noise_power = F_radar*0.001*10^(Lp/10)*Bw;
radar.Power = 10^(radar.SNR/10)*radar.noise_power*ones(radar.Tx,1);
radar.clutter_power = 10^(radar.CNR/10)*radar.noise_power;
%% Alternating optimization
radar.iota_max = 8;% algorithm 3
fdcomm.tu_max = 7; % algorithm 1 max number of iterations to execute the UL subgradient method
fdcomm.td_max = 7; % algorithm 2 % max number of iterations to execute the DL suggradient method
%% Subgradient method Algorithm 1 & 2
fdcomm.lambda_UL = ones(fdcomm.UL_num,radar.codelength);
fdcomm.lambda_DL = ones(radar.codelength,1);
fdcomm.mu_UL = ones(fdcomm.UL_num,radar.codelength);
fdcomm.mu_DL = ones(fdcomm.DL_num,radar.codelength);
%% comm rate
fdcomm.R_DL = systemCfg.Rates.DL;
fdcomm.R_UL = systemCfg.Rates.UL;
%% Radar comm symbols
radar_comm.n_Bm = systemCfg.Radar_Comm_symbols.n_Bm;
radar_comm.nu = systemCfg.Radar_Comm_symbols.n_u;
radar.CUT_Idx = systemCfg.Radar_Comm_symbols.n_t;
%% Covariance matrix initializations
cov.Bmr = zeros(radar.codelength,radar.codelength,radar.Rx);
cov.Btr = zeros(radar.codelength,radar.codelength,radar.Rx);
cov.radar2BS = zeros(fdcomm.BSRx,fdcomm.BSRx,radar.codelength);
cov.UL2radar = zeros(radar.codelength,radar.codelength,radar.Rx);

R_BJ = cell(fdcomm.DL_num,radar.codelength);
R_MUI_DL = cell(fdcomm.DL_num,radar.codelength);
R_IB = cell(fdcomm.UL_num,radar.codelength);
R_MUI_UL = cell(fdcomm.UL_num,radar.codelength);
R_BB = cell(radar.codelength,1);
R_ULDL = cell(fdcomm.DL_num,radar.codelength);
for k = 1:radar.codelength
    R_BB{k,1} = zeros(fdcomm.BSRx,fdcomm.BSTx);
    for jj = 1:fdcomm.DL_num
        R_BJ{jj,k} = zeros(fdcomm.DL_UE_Ant(jj),fdcomm.DL_UE_Ant(jj));
        R_MUI_DL{jj,k} = zeros(fdcomm.DL_UE_Ant(jj),fdcomm.DL_UE_Ant(jj));
        R_ULDL{jj,k} = zeros(fdcomm.DL_UE_Ant(jj),fdcomm.DL_UE_Ant(jj));
    end
    for ii = 1:fdcomm.UL_num
        R_IB{ii,k} = zeros(fdcomm.BSRx,fdcomm.BSRx);
        R_MUI_UL{ii,k} = zeros(fdcomm.BSRx,fdcomm.BSRx);
    end
end
cov.DL = R_BJ;
cov.MUI_DL = R_MUI_DL;
cov.MUI_UL = R_MUI_UL;
cov.UL = R_IB;
cov.B2B = R_BB;
cov.UL2DL = R_ULDL;
cov.in_DL = cov.DL;
cov.total_DL = cov.DL;
cov.in_UL = cov.UL;
cov.total_UL = cov.UL;
R_rJ = cell(fdcomm.DL_num,1);
for jj = 1:fdcomm.DL_num
    R_rJ{jj} = zeros(fdcomm.DL_UE_Ant(jj),radar.Tx);
end
cov.radar2DL = R_rJ;
R_Zr = zeros(radar.codelength,radar.codelength,radar.Rx);
for nr = 1:radar.Rx
    R_Zr(:,:,nr) = radar.noise_power*eye(radar.codelength);
end
cov.noise = R_Zr;
end