function [cov] = tsp_covmat_rev(fdcomm,radar,radar_comm,cov)

R_Bmr = cov.Bmr;
R_rB = cov.radar2BS;
R_U_Nr_total = cov.UL2radar;
R_BJ = cov.DL;
R_MUI_DL = cov.MUI_DL;
R_MUI_UL = cov.MUI_UL;
R_IB = cov.UL;
R_BB = cov.B2B;
R_ULDL = cov.UL2DL;
%% Covariances at radar Rxs
%Mr = radar.Tx; % Number of radar TX antennas
Nr = radar.Rx; % Number of radar RX antennas
K = radar.codelength; % The length of the radar code; or the number of PRIs
%nt = radar.CUT_Idx; % CUT index 
A = radar.codematrix;
R_rt_r = zeros(K,K,Nr);
Sigma_c_r = radar.Clutter_channel_cov; % Clutter channel covariance matrix
R_c_r = zeros(K,K,Nr);
Sigma_rt_Nr = radar.Sigma_rt_Nr;
S_tr_cell = cell(K,1);
for k = 1:K
    S_tr_cell{k,1} = A(k,:);
end
S_tr = blkdiag(S_tr_cell{:});
cov.S_tr = S_tr;
for nr = 1:Nr
    for m = 1:K
        a_m = A(m,:).';
        for l = 1:K
            a_l = A(l,:).';
            Sigma_rt_nr_m_l = Sigma_rt_Nr{m,l,nr};
            R_rt_r(m,l,nr)=abs(trace(a_m*a_l'*Sigma_rt_nr_m_l));
            R_c_r(m,l,nr) = abs(trace(a_m*a_l'*Sigma_c_r));
        end
    end
end
cov.rtr = R_rt_r;
R_tr = R_rt_r;
cov.cr = R_c_r;
cov.target2radar = cov.rtr;
%% FD comm related covariance matrices
if fdcomm.DL_num>0 % DL is enabled
    Mc = fdcomm.BSTx;% Number of BS TX antennas
    J = fdcomm.DL_num; % Number of DL UEs
    N_DL = fdcomm.DL_UE_Ant; % Number of the DL UE antennas
    D_DL = fdcomm.DLsymbols;
    P_dJ = fdcomm.DLprecoders; % cell(J,K*N)
    n_Bm = radar_comm.n_Bm;
    %% DL - radar
    S_Bm = zeros(Mc,K);
    for k = 1:K
        s_dj_k_nBm = 0;
        for jj = 1:J  
            P_Bj_k  = P_dJ{jj,k};
            s_dj_k_nBm = P_Bj_k*D_DL{jj}(:,n_Bm,k) + s_dj_k_nBm;
        end
        S_Bm(:,k) = s_dj_k_nBm;
    end
    R_Bmr = zeros(K,K,Nr);
    Sigma_Bm_Nr = radar_comm.Sigma_Bm_Nr;
    for nr = 1:Nr
        for m = 1:K
            S_Bm_m = S_Bm(:,m);
            for l = 1:K
                S_Bm_l = S_Bm(:,l);
                Sigma_Bm_nr_m_l = Sigma_Bm_Nr{m,l,nr};
                R_Bmr(m,l,nr)=abs(trace(S_Bm_m*S_Bm_l'*Sigma_Bm_nr_m_l));
            end
        end
    end
    cov.Bmr = R_Bmr;
    cov.S_Bm = S_Bm;
    %% BS-target-radar
    if radar_comm.isCollaborate
        S_Bt = zeros(Mc,K);
        S_t_cell = cell(K,1);
        %S_t_K = zeros(Mc+Mr,K);
        for k = 1:K
            s_dj_k_one = 0;
            a_k = A(k,:).';
            for jj = 1:J  
                P_Bj_k  = P_dJ{jj,k};
                s_dj_k_one  = P_Bj_k*D_DL{jj}(:,1,k) + s_dj_k_one;
            end
            S_Bt(:,k) = s_dj_k_one;
            s_t_k = [a_k;s_dj_k_one];
            %S_t_K(:,k) = s_t_k;
            S_t_cell{k,1} =s_t_k.';
        end
        % https://www.mathworks.com/matlabcentral/answers/46316-sparse-block-diagonal-matrix
        S_tr = blkdiag(S_t_cell{:}); 
        cov.S_tr = S_tr;
        R_Btr = zeros(K,K,Nr);
        R_tr = zeros(K,K,Nr);
        Sigma_Bt_Nr = radar_comm.Sigma_Bt_Nr;
        for nr = 1:Nr
            for m = 1:K
                S_Bt_m = S_Bt(:,m);
                for l = 1:K
                    S_Bt_l = S_Bt(:,l);
                    Sigma_Bt_nr_m_l = Sigma_Bt_Nr{m,l,nr};
                    R_Btr(m,l,nr)=abs(trace(S_Bt_m*S_Bt_l'*Sigma_Bt_nr_m_l));
                    R_tr(m,l,nr) = R_Btr(m,l,nr)+R_rt_r(m,l,nr);
                end
            end
        end
        cov.Btr = R_Btr;
        cov.target2radar = R_tr;
        cov.S_t_cell = S_t_cell;
        cov.S_Bt = S_Bt;
    end
    %% radar to DL UEs
    R_rJ    = cell(J,1);
    H_r_DL  = radar_comm.radar2DLchannnels;
    for jj = 1:J
        N_j = N_DL(jj);
        R_rj = zeros(N_j,N_j);
        H_rj = H_r_DL{jj};
        for k = 1:K
            ak = A(k,:).';
            R_rj_k = H_rj*(ak*ak')*H_rj';
            d = eye(size(R_rj_k), 'logical');
            R_rj_k(d) = real(diag(R_rj_k));
            R_rj(:,:,k) = R_rj_k;
        end
        R_rJ{jj,1} = R_rj;
    end
    cov.radar2DL = R_rJ;
    %% BS-DL UE
    H_DL = fdcomm.DLchannels;
    R_BJ = cell(J,K);   % All the DL covarianc matrices
    R_MUI_DL = cell(J,K);
    for k = 1 : K
        for ii = 1:J
            HBj = H_DL{ii}; %load the UL channel matrix
            PBj_k = P_dJ{ii,k};
            R_Bj_k = HBj*(PBj_k*PBj_k')*HBj'; 
            d = eye(size(R_Bj_k), 'logical');
            R_Bj_k(d) = real(diag(R_Bj_k));
            R_BJ{ii,k} = R_Bj_k;
        end
    end
    cov.DL= R_BJ;
    for k = 1 : K
        for jj = 1:J
            HBj = H_DL{jj}; %load the UL channel matrix
            jj_prime = [1:jj-1 jj+1:J];
            Nj = N_DL(jj);
            R_j_MUI = zeros(Nj,Nj);
            for jjj = 1:length(jj_prime)
                j_mui = jj_prime(jjj);
                Pj_mui = P_dJ{j_mui,k};
                R_jjj_MUI = HBj*(Pj_mui*Pj_mui')*HBj';
                d = eye(size(R_jjj_MUI), 'logical');
                R_jjj_MUI(d) = real(diag(R_jjj_MUI));
                R_j_MUI = R_jjj_MUI +R_j_MUI;
            end
            R_MUI_DL{jj,k} = R_j_MUI;
        end
    end
    cov.MUI_DL = R_MUI_DL;
end
%% UL related
if fdcomm.UL_num>0
    I = fdcomm.UL_num; % Number of UL UEs
    P_uI = fdcomm.ULprecoders; % cell(I,k*N)
    D_UL = fdcomm.ULsymbols;
    Nc = fdcomm.BSRx;% Number of BS RX antennas
    % S_U_r
    nu = radar_comm.nu;
    S_Ur = cell(K,I,Nr);
    for nr = 1:Nr
        for ii =1:I
            for k = 1:K
                P_ui_k = P_uI{ii,k};
                d_ui_k = D_UL{ii}(:,nu,k);
                S_Ur{k,ii,nr} = P_ui_k*d_ui_k;
            end
        end
    end
    Sigma_U_Nr = radar_comm.Sigma_U_Nr;
    R_U_Nr = cell(I,Nr);
    R_U_Nr_total = zeros(K,K,Nr);
    for nr = 1:Nr
        R_U_nr_temp = 0;
        for ii = 1:I
            R_i_nr = zeros(K,K);
            for m = 1:K
                s_inr_m = S_Ur{m,ii,nr};
               for l = 1:K
                   s_inr_l = S_Ur{l,ii,nr};
                   R_i_nr(m,l) = abs(trace(s_inr_m*s_inr_l'*Sigma_U_Nr{m,l,ii,nr}));
               end
            end
            R_U_Nr{ii,nr} = R_i_nr;
            R_U_nr_temp = R_U_nr_temp + R_i_nr;
       end
       R_U_Nr_total(:,:,nr) = R_U_nr_temp;
    end
    cov.S_Ur = S_Ur;
    cov.U_r = R_U_Nr;
    cov.UL2radar = R_U_Nr_total;
    %% radar to BS
    R_rB = zeros(Nc,Nc,K);
    H_r_BS = radar_comm.radar2BSchannels;
    for k = 1:K 
        ak = A(k,:).';
        R_rB_k = H_r_BS*(ak*ak')*H_r_BS';
        d = eye(size(R_rB_k), 'logical');
        R_rB_k(d) = abs(diag(R_rB_k));
        R_rB(:,:,k) = R_rB_k;
    end
    cov.radar2BS = R_rB;
    %% UL UE - BS
    H_UL = fdcomm.ULchannels;
    R_IB = cell(I,K);   % All the UL covariance matrices
    R_MUI_UL = cell(I,K);
    R_sum_UL = cell(K,1);
    for k = 1:K
        R_sum = zeros(Nc,Nc);
        for ii = 1 : I
            H_iB = H_UL{ii}; %load the UL channel matrix
            PiB_k = P_uI{ii,k};
            R_iB_k = H_iB*(PiB_k*PiB_k')*H_iB';
            d = eye(size(R_iB_k), 'logical');
            R_iB_k(d) = abs(diag(R_iB_k));
            R_IB{ii,k} = R_iB_k;
            R_sum = R_IB{ii,k}+R_sum;
        end
        R_sum_UL{k} = R_sum;
    end
    for k = 1:K
        R_sum_k = R_sum_UL{k};
        for ii = 1:I
            R_MUI_UL_ii_k = R_sum_k - R_IB{ii,k};
            d = eye(size(R_MUI_UL_ii_k), 'logical');
            R_MUI_UL_ii_k(d) = abs(diag(R_MUI_UL_ii_k));
            R_MUI_UL{ii,k} = R_MUI_UL_ii_k;
        end
    end
    cov.MUI_UL = R_MUI_UL;
    cov.UL = R_IB;
    cov.sum_UL = R_sum_UL;
end
%% FD enabled
if fdcomm.UL_num>0 && fdcomm.DL_num > 0
    H_UL_DL = fdcomm.ULDLchannels;
    H_BB = fdcomm.BBchannel;
    R_BB = cell(K,1);
    % BS - BS 
    for k = 1:K
        R_BB_k = H_BB*(S_Bm(:,k)*S_Bm(:,k)')*H_BB';
        d = eye(size(R_BB_k), 'logical');
        R_BB_k(d) = abs(diag(R_BB_k));
        R_BB{k,1} = R_BB_k; 
    end
    cov.B2B = R_BB;
    % UL to DL
    R_ULDL = cell(J,K);
    for jj = 1:J
        for kk = 1:K
            R_ULDL_temp = 0; 
            for ii = 1:I
                Hij = H_UL_DL{ii,jj}; %load the UL channel matrix
                PiB_k = P_uI{ii,k};
                R_ij = Hij*(PiB_k*PiB_k')*Hij';
                d = eye(size(R_ij), 'logical');
                R_ij(d) = abs(diag(R_ij));
                R_ULDL_temp = R_ij +R_ULDL_temp;
            end
            R_ULDL{jj,kk} = R_ULDL_temp;
        end
    end
    cov.UL2DL = R_ULDL;
end

%% Radar Total Cov
% radar interference matrices
R_Zr = cov.noise;
R_in_r = zeros(K,K,Nr);
R_r = zeros(K,K,Nr);
for nr = 1:Nr
    R_in_r(:,:,nr) = R_Zr(:,:,nr)+ R_U_Nr_total(:,:,nr)+ R_Bmr(:,:,nr)+ R_c_r(:,:,nr);
    R_r(:,:,nr) = R_in_r(:,:,nr)+R_tr(:,:,nr);
end
cov.inr = R_in_r;
cov.total_r = R_r;

%% UL IN and total covs
if fdcomm.UL_num>0
    R_in_UL = cell(I,K);
    R_total_UL = cell(I,K);
    R_ZB = fdcomm.BS_noise_power*eye(Nc);
    for ii = 1:I
        for k = 1:K
            R_MUI_i_k = R_MUI_UL{ii,k};
            R_BB_k = R_BB{k,1};
            R_rB_k = R_rB(:,:,k);
            R_iB_k = R_IB{ii,k};
            R_in_UL{ii,k} = R_MUI_i_k+ R_BB_k + R_rB_k + R_ZB;
            R_total_UL{ii,k} = R_in_UL{ii,k} + R_iB_k;
        end
    end
    cov.in_UL = R_in_UL;
    cov.total_UL = R_total_UL;
end

%% DL total Covs
if fdcomm.DL_num>0
    R_in_DL = cell(J,K);
    R_total_DL = cell(J,K);
    for jj = 1:J
        R_rj = R_rJ{jj,1};
        N_j = N_DL(jj);
        for k = 1:K
            R_MUI_j_k = R_MUI_DL{jj,k};
            R_UL_j_k  = R_ULDL{jj,k};
            R_rj_k = R_rj(:,:,k);
            R_in_DL{jj,k} = R_MUI_j_k + R_UL_j_k + R_rj_k + fdcomm.UE_noise_power*eye(N_j);
            R_dj_k = R_BJ{jj,k};
            R_total_DL{jj,k} = R_dj_k + R_in_DL{jj,k};
        end
    end
    cov.in_DL = R_in_DL;
    cov.total_DL = R_total_DL;
end
end

