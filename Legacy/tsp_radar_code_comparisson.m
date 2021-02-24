%%% Convergence analysis

%% 
SNR_Btr = -20:5:10;
file_name = ['uniformweight_BTR_', num2str(snr_Btr),'dB_co.mat'];
%% Initialization radar code are uncoded using random initialization
[fdcomm, radar, radar_comm, cov] = tsp_ini_random(SNR,radar,fdcomm);
radar_op = radar;
fdcomm_op = fdcomm; 
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_max = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
while ell <= ell_max
    disp(ell);
    %% WMMSE
    iota = 1;
    while iota<= iota_max
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            %disp(['ii=',num2str(ii)]);
            [fdcomm] = tsp_UL_napla(k, ii, fdcomm, radar, radar_comm); %% modify the 
            [fdcomm,radar,cov] = tsp_UL_subgradient(fdcomm, radar_comm, radar, cov, ii,k);
        end
        %% DL precoder
        for jj = 1 : J
            %disp(['jj=',num2str(ii)])
            [fdcomm] = tsp_DL_napla(k, jj, fdcomm, radar, radar_comm, cov);
            [fdcomm,radar,cov] = tsp_DL_subgradient(fdcomm, radar_comm, radar, cov, jj ,k);
        end
    end
    iota = iota + 1;
    end
    % Update covariance matrices
    cov = tsp_covmat(fdcomm,radar,radar_comm);
    % Update linear receivers
    fdcomm = tsp_Comm_MMSE(fdcomm, radar, cov);
    radar = tsp_radar_MMSE(radar, cov);
    I_total(ell) = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_max = I_total(ell);
        radar_op = radar;
        fdcomm_op = fdcomm;
    else
        I_total_op(ell) = I_max;
    end
    ell = ell+1;
end
save('uniformcoded.mat');
%% Initialization radar code are uncoded using random initialization
[fdcomm, radar, radar_comm, cov] = tsp_ini_random(SNR,radar,fdcomm);
radar.codematrix = 
radar_op = radar;
fdcomm_op = fdcomm; 
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_max = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
while ell <= ell_max
    disp(ell);
    %% WMMSE
    iota = 1;
    while iota<= iota_max
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            %disp(['ii=',num2str(ii)]);
            [fdcomm] = tsp_UL_napla(k, ii, fdcomm, radar, radar_comm); %% modify the 
            [fdcomm,radar,cov] = tsp_UL_subgradient(fdcomm, radar_comm, radar, cov, ii,k);
        end
        %% DL precoder
        for jj = 1 : J
            %disp(['jj=',num2str(ii)])
            [fdcomm] = tsp_DL_napla(k, jj, fdcomm, radar, radar_comm, cov);
            [fdcomm,radar,cov] = tsp_DL_subgradient(fdcomm, radar_comm, radar, cov, jj ,k);
        end
    end
    iota = iota + 1;
    end
    % Update covariance matrices
    cov = tsp_covmat(fdcomm,radar,radar_comm);
    % Update linear receivers
    fdcomm = tsp_Comm_MMSE(fdcomm, radar, cov);
    radar = tsp_radar_MMSE(radar, cov);
    I_total(ell) = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_max = I_total(ell);
        radar_op = radar;
        fdcomm_op = fdcomm;
    else
        I_total_op(ell) = I_max;
    end
    ell = ell+1;
end

