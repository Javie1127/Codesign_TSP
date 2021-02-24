function [fdcomm, radar] = tsp_ini_mmse(radar,fdcomm,cov)
%%%%-----------------------------------
%% Initializing the MMSE matrices
K = radar.codelength;
fdcomm = tsp_Comm_MMSE(fdcomm,radar,cov);
radar = tsp_radar_MMSE(radar,cov);
%% Initializing the performance measures
radar = Xi_radar(radar);
for k = 1:K
    fdcomm = Xi_comm_k(fdcomm,k);
end
end
