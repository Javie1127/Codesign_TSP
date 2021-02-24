%%% Convergence analysis for equal weights

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
% SNR.rtr = randi([-5,5],Mr,Nr);
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
% Clutter 
SNR.CNR = ones(Nr,1);
ell_max = 20; % algorithm 5
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);

%% simulation parameters
[fdcomm_para_ew, radar_para_ew, radar_comm_para] = tsp_parameters(SNR,radar,fdcomm);
%% random initialization
M = 15;
I_total_ew_ri = zeros(ell_max,M);
parfor m = 1:M
    %% Random Initialization
    [fdcomm_m_ew_ri, radar_m_ew_ri, cov_m_ew_ri] = tsp_ini_random(radar_para_ew,fdcomm_para_ew,radar_comm_para);
    %% Alternating optimization
    % radar_m_op_ew_ri = radar_m_ini_ew_ri;
    % fdcomm_m_op_ew_ri = fdcomm_m_ini_ew_ri;
    I_total_m_ew_ri = zeros(ell_max,1);
    I_total_m_op_ew_ri = zeros(ell_max,1);
    I_max_m_ew_ri = sum(fdcomm_m_ew_ri.MI_UL(:))+sum(fdcomm_m_ew_ri.MI_DL(:))+ sum(radar_m_ew_ri.MI_radar);
    ell = 1;
    while ell <= ell_max
        [fdcomm_m_ew_ri, radar_m_ew_ri_ast] =...
            tsp_WMMSE_algorithm(fdcomm_m_ew_ri, radar_m_ew_ri, radar_comm_para, cov_m_ew_ri);
        % Calculate A^star
        radar_m_ew_ri = tsp_Nearest_PAR(radar_m_ew_ri_ast);
        % Update covariance matrices
        cov_m_ew_ri = tsp_covmat(fdcomm_m_ew_ri,radar_m_ew_ri,radar_comm_para);
        % Update linear receivers
        fdcomm_m_ew_ri = tsp_Comm_MMSE(fdcomm_m_ew_ri, radar_m_ew_ri, cov_m_ew_ri);
        radar_m_ew_ri = tsp_radar_MMSE(radar_m_ew_ri, cov_m_ew_ri);
        I_total_m_ew_ri(ell) = sum(fdcomm_m_ew_ri.MI_UL(:))+sum(fdcomm_m_ew_ri.MI_DL(:))...
            + sum(radar_m_ew_ri.MI_radar);
        if I_total_m_ew_ri(ell) > I_max_m_ew_ri
            I_total_m_op_ew_ri(ell) = I_total_m_ew_ri(ell);
            I_max_m_ew_ri = I_total_m_ew_ri(ell);
            %radar_m_op_ew_ri = radar_m_ew_ri;
            %fdcomm_m_op_ew_ri = fdcomm_m_ew_ri;
        else
            I_total_m_op_ew_ri(ell) = I_max_m_ew_ri;
        end
        ell = ell+1;
    end
    %fdcomm_m_op_ew_ri.I_total = I_total_m_ew_ri;
    %fdcomm_m_op_ew_ri.I_total_op = I_total_m_op_ew_ri;
    %cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
    I_total_ew_ri(:,m) = I_total_m_op_ew_ri;
end
I_total_avg_ew_ri = mean(I_total_ew_ri,2);
%% normal initialization
[fdcomm_ini_ew_ni, radar_ini_ew_ni, cov_ew_ni] = tsp_ini_normal(radar_para_ew,fdcomm_para_ew,radar_comm_para);
%% Alternating optimization

% fdcomm_op_ew_ni = fdcomm_ini_ew_ni;
% radar_op_ew_ni = radar_ini_ew_ni;
fdcomm_ew_ni = fdcomm_ini_ew_ni;
radar_ew_ni = radar_ini_ew_ni;
I_total_ew_ni = zeros(ell_max,1);
I_total_op_ew_ni = zeros(ell_max,1);
I_max_ew_ni = sum(fdcomm_ew_ni.MI_UL(:))+sum(fdcomm_ew_ni.MI_DL(:))+ sum(radar_ew_ni.MI_radar);
ell = 1;
while ell <= ell_max
    [fdcomm_ew_ni, radar_ew_ast] = ...
        tsp_WMMSE_algorithm(fdcomm_ew_ni, radar_ew_ni, radar_comm_para, cov_ew_ni);
    % Calculate A^star
    radar_ew_ni = tsp_Nearest_PAR(radar_ew_ast);
    % Update covariance matrices
    cov_ew_ni = tsp_covmat(fdcomm_ew_ni,radar_ew_ni,radar_comm_para);
    % Update linear receivers
    fdcomm_ew_ni = tsp_Comm_MMSE(fdcomm_ew_ni, radar_ew_ni, cov_ew_ni);
    radar_ew_ni = tsp_radar_MMSE(radar_ew_ni, cov_ew_ni);
    I_total_ew_ni(ell) = sum(fdcomm_ew_ni.MI_UL(:))+sum(fdcomm_ew_ni.MI_DL(:))+ sum(radar_ew_ni.MI_radar);
    if I_total_ew_ni(ell) > I_max_ew_ni
        I_total_op_ew_ni(ell) = I_total_ew_ni(ell);
        I_max_ew_ni = I_total_ew_ni(ell);
        % radar_op = radar;
    else
        I_total_op_ew_ni(ell) = I_max_ew_ni;
    end
    ell = ell+1;
end
% fdcomm_op_ew_ni.I_total = I_total;
% fdcomm_op_ew_ni.I_total_op = I_total_op;
%cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);


fdcomm_para_uew = fdcomm_para_ew;
radar_para_uew = radar_para_ew;
%% priority unequal weight
fdcomm_para_uew.alpha_UL = 0.1*ones(I,1);
fdcomm_para_uew.alpha_DL = 0.05*ones(J,1);
radar_para_uew.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% random initialization
I_total_uew_ri = zeros(ell_max,M);
parfor m = 1:M
    %% Random Initialization
    [fdcomm_m_uew_ri, radar_m_uew_ri, cov_m_uew_ri] = tsp_ini_random(radar_para_uew,fdcomm_para_uew,radar_comm_para);
    %% Alternating optimization
    %radar_op = radar_ini;
%     fdcomm_op_uew_ri = fdcomm_ini;
    I_total_m_uew_ri = zeros(ell_max,1);
    I_total_op_m_uew_ri = zeros(ell_max,1);
    I_max_m_op_uew_ri = sum(fdcomm_m_uew_ri.MI_UL(:))+sum(fdcomm_m_uew_ri.MI_DL(:))+ sum(radar_m_uew_ri.MI_radar);
    ell = 1;
    while ell <= ell_max
        [fdcomm_m_uew_ri, radar_m_uew_ri_ast] =...
            tsp_WMMSE_algorithm(fdcomm_m_uew_ri, radar_m_uew_ri, radar_comm_para, cov_m_uew_ri);
        % Calculate A^star
        radar_m_uew_ri = tsp_Nearest_PAR(radar_m_uew_ri_ast);
        % Update covariance matrices
        cov_m_uew_ri = tsp_covmat(fdcomm_m_uew_ri,radar_m_uew_ri,radar_comm_para);
        % Update linear receivers
        fdcomm_m_uew_ri = tsp_Comm_MMSE(fdcomm_m_uew_ri, radar_m_uew_ri, cov_m_uew_ri);
        radar_m_uew_ri = tsp_radar_MMSE(radar_m_uew_ri, cov_m_uew_ri);
        I_total_m_uew_ri(ell) = sum(fdcomm_m_uew_ri.MI_UL(:))+sum(fdcomm_m_uew_ri.MI_DL(:))+ sum(radar_m_uew_ri.MI_radar);
        if I_total_m_uew_ri(ell) > I_max_m_op_uew_ri
            I_total_op_m_uew_ri(ell) = I_total_m_uew_ri(ell);
            I_max_m_op_uew_ri = I_total_m_uew_ri(ell);
            % radar_op = radar;
%             fdcomm_op_uew_ri = fdcomm;
        else
            I_total_op_m_uew_ri(ell) = I_max_m_op_uew_ri;
        end
        ell = ell+1;
    end
%     fdcomm_op_uew_ri.I_total = I_total_m_uew_ri;
%     fdcomm_op_uew_ri.I_total_op = I_total_op;
    % cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
    I_total_uew_ri(:,m) = I_total_op_m_uew_ri;
end
I_total_avg_uew_ri = mean(I_total_uew_ri,2);

%% normal initialization
[fdcomm_uew_ni, radar_uew_ni, cov_uew_ni] = tsp_ini_normal(radar_para_uew,fdcomm_para_uew,radar_comm_para);
%% Alternating optimization
% radar_op = radar_ini;
% fdcomm_op_uew_ni = fdcomm_ini;
% radar_uew_ni = radar_ini;
% fdcomm = fdcomm_ini;
I_total_uew_ni = zeros(ell_max,1);
I_total_op_uew_ni = zeros(ell_max,1);
I_max_uew_ni = sum(fdcomm_uew_ni.MI_UL(:))+sum(fdcomm_uew_ni.MI_DL(:))+ sum(radar_uew_ni.MI_radar);
ell = 1;
while ell <= ell_max
    [fdcomm_uew_ni, radar_uew_ni] = ...
        tsp_WMMSE_algorithm(fdcomm_uew_ni, radar_uew_ni, radar_comm_para, cov_uew_ni);
    % Calculate A^star
    radar_uew_ni = tsp_Nearest_PAR(radar_uew_ni);
    % Update covariance matrices
    cov_uew_ni = tsp_covmat(fdcomm_uew_ni,radar_uew_ni,radar_comm_para);
    % Update linear receivers
    fdcomm_uew_ni = tsp_Comm_MMSE(fdcomm_uew_ni, radar_uew_ni, cov_uew_ni);
    radar_uew_ni = tsp_radar_MMSE(radar_uew_ni, cov_uew_ni);
    I_total_uew_ni(ell) = sum(fdcomm_uew_ni.MI_UL(:))+sum(fdcomm_uew_ni.MI_DL(:))+ sum(radar_uew_ni.MI_radar);
    if I_total_uew_ni(ell) > I_max_uew_ni
        I_total_op_uew_ni(ell) = I_total_uew_ni(ell);
        I_max_uew_ni = I_total_uew_ni(ell);
    else
        I_total_op_uew_ni(ell) = I_max_uew_ni;
    end
    ell = ell+1;
end
 
