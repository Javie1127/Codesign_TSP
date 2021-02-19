function [cov] = tsp_covmat(fdcomm,radar,radar_comm)

%% radar only covariances
Mr = radar.TX; % Number of radar TX antennas
Nr = radar.RX; % Number of radar RX antennas
K = radar.codelength; % The length of the radar code; or the number of PRIs
n = radar.CUT_Idx; % CUT index
%H_rtr                       = radar.channel;
Qr = radar.doppler; %Mr*K*Nr radar temporal steering vector
R_C = radar.cluttercov;
eta_rtr = radar.channelgain;
A = radar.codematrix;
sigma_0 = radar.noisepower;
%% FD comm related covariance matrices
Mc = fdcomm.BSTX;% Number of BS TX antennas
Nc = fdcomm.BSRX;% Number of BS RX antennas
I = fdcomm.UL_num; % Number of UL UEs
J = fdcomm.DL_num; % Number of DL UEs
N_DL = fdcomm.DL_UE_Ant; % Number of the DL UE antennas
P_dJ = fdcomm.DLprecoders; % cell(J,K*N)
P_uI = fdcomm.ULprecoders; % cell(I,k*N)
H_UL_DL = fdcomm.ULDLchannels;
D_DL = fdcomm.DLsymbols;
D_UL = fdcomm.ULsymbols;
H_DL = fdcomm.DLchannels;
H_UL = fdcomm.ULchannels;
H_BB = fdcomm.BBchannel;
% eta_BB = fdcomm.BBchannelgain;
R_IB = cell(I,K);   % All the UL covariance matrices
R_BJ = cell(J,K);   % All the DL covarianc matrices
R_MUI_DL = cell(J,K);
R_MUI_UL = cell(I,K);
% R_sum_DL = cell(K,1);
R_sum_UL = cell(K,1);
R_BB = cell(K,1);
%% radar-comm 
% H_Bmr = radar_comm.Bmrchannels;
% H_Btr = radar_comm.Btrchannels;
eta_Bmr = radar_comm.Bmrchannelgains;
eta_Btr = radar_comm.Btrchannelgains;
eta_UL_r = radar_comm.UL2rchannelgains;
% H_UL_r = radar_comm.UL2rchannles;
f_Bt_Nr = radar_comm.f_Bt_Nr;
f_Bm_Nr = radar_comm.f_Bm_Nr;
% Q_Btr = radar_comm.Btrdoppler;
% Q_Bmr = radar_comm.Bmrdoppler;
Q_Ir = radar_comm.ULrdoppler;
%% Covariance matrix at the nr radar RX 
% DL signals
S_dJ_one = zeros(Mc,K);
S_dJ_n = zeros(Mc,K);
for k = 1:K
    s_dj_k_one = 0;
    s_dj_k_n = 0;
    for jj = 1:J  
        P_Bj_k  = P_dJ{jj,k};
        s_dj_k_one  = P_Bj_k*D_DL{jj}(:,1,k) + s_dj_k_one;
        s_dj_k_n = P_Bj_k*D_DL{jj}(:,n+1,k) + s_dj_k_n;
    end
    S_dJ_one(:,k) = s_dj_k_one;
    S_dJ_n(:,k) = s_dj_k_n;
end
% BS - BS 
for k = 1:K
    R_BB{k,1} = H_BB*(S_dJ_n(:,k)*S_dJ_n(:,k)')*H_BB';
end

% covariance matrix R_t_r 
% BS-multipath/target-radar
R_Bmr = zeros(K,K,Nr);
for nr = 1:Nr
    f_Bmnr = f_Bm_Nr(nr);
    eta_Bmnr = eta_Bmr(nr);
    for m = 1:K
        for l = 1:K
            R_Bmr(m,l,nr) = eta_Bmnr*exp(1i*2*pi*(k-1)*f_Bmnr)*S_dJ_n(:,m).'*S_dJ_n(:,l);
        end
    end
   R_Bmr(:,:,nr) = R_Bmr(:,:,nr); 
end
S_rtr = zeros(K,Mr,Nr); % signal matrix of the radar 
S_Btr = zeros(K,Mc,Nr);
M = Mc+Mr;
S_tr = zeros(K,M,Nr);
for nr = 1:Nr
    f_Btnr = f_Bt_Nr(nr);
    for k = 1:K
        q_rnr_k = Qr(:,k,nr);
        Q_rnr_k = diag(q_rnr_k);
        s_Bt_k = exp(1i*2*pi*(k-1)*f_Btnr)*S_dJ_one(:,k);
        a_k = A(k,:).';
        s_rtr_k = Q_rnr_k*a_k;
        s_tr_k = [s_rtr_k.',s_Bt_k.'].';
        S_rtr(k,:,nr) = s_rtr_k.';
        S_Btr(k,:,nr) = s_Bt_k.';
        S_tr(k,:,nr) = s_tr_k.';
    end
end
R_tr = zeros(K,K,Nr); %Covariance for the radar
for nr = 1:Nr
    eta_rtnr = eta_rtr(nr);
    eta_Btnr = eta_Btr(nr);
    Sigma_tnr = radar.Sigma(:,:,nr);
    for m = 1:K 
        s_rt_nr_m = S_rtr(m,:,nr).';
        s_Bt_nr_m = S_Btr(m,:,nr).';
        s_tnr_m = S_tr(m,:,nr).';
        for l = 1:K
            s_rt_nr_l = S_rtr(l,:,nr).';
            s_Bt_nr_l = S_Btr(l,:,nr).';
            R_tr(m,l,nr) = eta_rtnr*s_rt_nr_l'*s_rt_nr_m+eta_Btnr*s_Bt_nr_l'*s_Bt_nr_m;
            s_tnr_l = S_tr(l,:,nr).';
            R_tr(m,l,nr) = trace(s_tnr_m*s_tnr_l'*Sigma_tnr);
        end
    end
    R_tr(:,:,nr) = R_tr(:,:,nr);
end
cov.S_tr = S_tr;
cov.S_rtr = S_rtr;
cov.S_Btr = S_Btr;
cov.target2radar = R_tr;
% UL UEs to radar RX
%H_UL_r = radar_comm.UL2radarchannels;
R_Ur = zeros(K,K,Nr);
for nr = 1:Nr
    for m= 1:K
        for l = 1:K
            R_U_nr_temp = 0;
            for ii = 1:I
                P_ui_m = P_uI{ii,m};
                d_ui_m = D_UL{ii}(:,1+n,m);
                s_inr_m = P_ui_m*d_ui_m;
                f_i_nr = log(Q_Ir(2,nr,ii))/(2*pi*1i);
                P_ui_l = P_uI{ii,l};
                d_ui_l = D_UL{ii}(:,1+n,l);
                s_inr_l = P_ui_l*d_ui_l;
                R_U_nr_temp = eta_UL_r(ii,nr)*exp(1i*2*pi*(m-l)*f_i_nr)*s_inr_l'*s_inr_m + R_U_nr_temp;
            end
            R_Ur(m,l,nr) = R_U_nr_temp;
        end
    end
    R_Ur(:,:,nr) = R_Ur(:,:,nr);
end
cov.UL2radar = R_Ur;
% radar interference matrices
R_Zr = zeros(K,K,Nr);
R_in_r = zeros(K,K,Nr);
sigma_n = radar.noisepower;
R_r = zeros(K,K,Nr);
for nr = 1:Nr
    R_Zr(:,:,nr) = sigma_n*eye(K,K);
    R_in_r(:,:,nr) = R_Zr(:,:,nr)+ R_Ur(:,:,nr)+ R_Bmr(:,:,nr)+ R_C(:,:,nr);
    R_r(:,:,nr) = R_in_r(:,:,nr)+R_tr(:,:,nr);
end
cov.noise = R_Zr;
cov.inr = R_in_r;
cov.total_r = R_r;
%% UL covs
for k = 1:K
    R_sum = zeros(Nc,Nc);
    for ii = 1 : I
        H_iB = H_UL{ii}; %load the UL channel matrix
        PiB_k = P_uI{ii,k};
        R_IB{ii,k} = H_iB*(PiB_k*PiB_k')*H_iB';
        R_sum = R_IB{ii,k}+R_sum;
    end
    R_sum_UL{k} = R_sum; 
end
for k = 1:K
    R_sum_k = R_sum_UL{k};
    for ii = 1:I
        R_MUI_UL_ii_k = R_sum_k - R_IB{ii,k};
        R_MUI_UL{ii,k} = R_MUI_UL_ii_k;
    end
end
cov.MUI_UL = R_MUI_UL;
% Self interference

% radar to BS
R_rB = zeros(Nc,Nc,K);
H_r_BS = radar_comm.radar2BSchannels;
for k = 1:K 
    ak = A(k,:).';
    R_rB(:,:,k) = H_r_BS*(ak*ak')*H_r_BS';
end
cov.radar2BS = R_rB;
% UL interference matrices
R_in_UL = cell(I,K);
R_total_UL = cell(I,K);
R_ZB = sigma_0*eye(Nc,Nc);
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
%% DL covs
for k = 1 : K
    for ii = 1:J
        HBj = H_DL{ii}; %load the UL channel matrix
        PBj_k = P_dJ{ii,k};
        R_BJ{ii,k} = HBj*(PBj_k*PBj_k')*HBj'; 
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
            R_j_MUI = HBj*(Pj_mui*Pj_mui')*HBj'+R_j_MUI;
        end
        R_MUI_DL{jj,k} = R_j_MUI;
    end
end
cov.MUI_DL = R_MUI_DL;
% for k = 1 : K
%     for jj = 1:J
%         R_bb = zeros(Nc,Nc);
%         PBj_k = P_dJ{jj,k};
%         R_bb = H_BB*(PBj_k*PBj_k')*H_BB'+R_bb;
%     end
%     R_BB{k,1} = nearestSPD(R_bb);
% end
cov.B2B = R_BB;
% UL to DL
R_ULDL = cell(J,I);
for jj = 1:J
    for kk = 1:K
        R_ULDL_temp = 0; 
        for ii = 1:I
            Hij = H_UL_DL{ii,jj}; %load the UL channel matrix
            PiB_k = P_uI{ii,k};
            R_ULDL_temp = Hij*(PiB_k*PiB_k')*Hij'+R_ULDL_temp;
        end
        R_ULDL{jj,kk} = R_ULDL_temp;
    end
end
cov.UL2DL = R_ULDL;
% radar to DL UEs
R_rJ    = cell(J,1);
H_r_DL  = radar_comm.radar2DLchannnels;
for jj = 1:J
    N_j = N_DL(jj);
    R_rj = zeros(N_j,N_j);
    H_rj = H_r_DL{jj};
    for k = 1:K
        ak = A(k,:).';
        R_rj(:,:,k) = H_rj*(ak*ak')*H_rj';
    end
    R_rJ{jj,1} = R_rj;
end
cov.radar2DL = R_rJ;
% DL interference matrices
R_in_DL = cell(J,K);
R_total_DL = cell(J,K);
for jj = 1:J
    R_rj = R_rJ{jj,1};
    N_j = N_DL(jj);
    for k = 1:K
        R_MUI_j_k = R_MUI_DL{jj,k};
        R_UL_j_k  = R_ULDL{jj,k};
        R_rj_k = R_rj(:,:,k);
        R_in_DL{jj,k} = R_MUI_j_k + R_UL_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
        R_dj_k = R_BJ{jj,k};
        R_total_DL{jj,k} = R_dj_k + R_MUI_j_k + R_UL_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
    end
end
cov.in_DL = R_in_DL;
cov.total_DL = R_total_DL;
end

