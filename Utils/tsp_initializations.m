function [fdcomm, radar,radar_comm] = tsp_initializations(systemCfg)
%% Radar-only parameters
radar.Tx =systemCfg.Num_Radar_Txs;
radar.RX =systemCfg.Num_Radar_Rxs;
radar.codelength = systemCfg.SLS_Top_Cfg.Num_PRI;
radar.num_range_cell = systemCfg.SLS_Top_Cfg.Num_Range_Cell;
radar.gamma_r = 2*ones(Mr,1); % PAR level
radar.K_factor = systemCfg.Channel_Modeling.K_Rician_radar;
radar.Rician_direct = systemCfg.Channel_Modeling.Rician_direct_radar;
%% FD-only parameters
if systemCfg.Cases.FD.isEnabled
    fdcomm.UL_num = systemCfg.Cases.FD.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.FD.Num_DL_UE;
    fdcomm.UL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas;
    fdcomm.DL_UE_Ant = fdcomm.UL_UE_Ant;
    fdcomm.ULstream_num = fdcomm.UL_UE_Ant;
    fdcomm.DLstream_num = fdcomm.DL_UE_Ant;
    radar_comm.isCollaborate = systemCfg.Cases.FD.isCollaborationEnabled;
    fdcomm.theta_BT = systemCfg.Channel_Modeling.theta_Bt;
elseif systemCfg.Cases.DL.isEnabled
    fdcomm.UL_num = systemCfg.Cases.DL.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.DL.Num_DL_UE;
    fdcomm.DL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas;
    fdcomm.UL_UE_Ant = 0;
    fdcomm.ULstream_num = fdcomm.UL_UE_Ant;
    fdcomm.DLstream_num = fdcomm.DL_UE_Ant;
    radar_comm.isCollaborate = systemCfg.Cases.DL.isCollaborationEnabled;
elseif systemCfg.Cases.UL.isEnabled
    fdcomm.UL_num = systemCfg.Cases.UL.Num_UL_UE;
    fdcomm.DL_num = systemCfg.Cases.UL.Num_DL_UE;
    fdcomm.DL_UE_Ant = 0;
    fdcomm.UL_UE_Ant = systemCfg.Antennas.Num_UE_Antennas;
    fdcomm.ULstream_num = fdcomm.UL_UE_Ant;
    fdcomm.DLstream_num = fdcomm.DL_UE_Ant;
    radar_comm.isCollaborate = 0;
end
fdcomm.pathloss = systemCfg.Channel_Modeling.Path_Loss;
fdcomm.K_factor = systemCfg.Channel_Modeling.K_Rician_SI;
fdcomm.UL_pathloss = systemCfg.Powers.SNR_UL;
% noise power
Lp = systemCfg.SLS_Top_Cfg.Noise_density;
Bw = systemCfg.SLS_Top_Cfg.System_Bandwidth;
NF_BS = systemCfg.SLS_Top_Cfg.BS_Noise_Figure;
F_BS = 10^((NF_BS)/10);
fdcomm.BS_noise_power = F_BS*0.001*10^(Lp/10)*Bw;
NF_UE = systemCfg.SLS_Top_Cfg.UE_Noise_Figure;
F_UE = 10^((NF_UE)/10);
fdcomm.BS_noise_power = F_UE*0.001*10^(Lp/10)*Bw;
NF_radar = systemCfg.SLS_Top_Cfg.Radar_Noise_Figure;
F_radar = 10^((NF_radar)/10);
radar.noise_power = F_radar*0.001*10^(Lp/10)*Bw;
end