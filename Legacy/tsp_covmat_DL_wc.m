function [cov] = tsp_covmat_DL_wc(DLcomm,radar,radar_comm)

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
%% FD comm related covariance matrices
Mc = DLcomm.BSTX;% Number of BS TX antennas
J = DLcomm.DL_num; % Number of DL UEs
N_DL = DLcomm.DL_UE_Ant; % Number of the DL UE antennas
P_dJ = DLcomm.DLprecoders; % cell(J,K*N)
D_DL = DLcomm.DLsymbols;
H_DL = DLcomm.DLchannels;
% eta_BB = fdcomm.BBchannelgain;
R_BJ = cell(J,K);   % All the DL covarianc matrices
R_MUI_DL = cell(J,K);
% R_sum_DL = cell(K,1);
%% radar-comm 
% H_Bmr = radar_comm.Bmrchannels;
% H_Btr = radar_comm.Btrchannels;
eta_Bmr = radar_comm.Bmrchannelgains;
% H_UL_r = radar_comm.UL2rchannles;
f_Bm_Nr = radar_comm.f_Bm_Nr;
% Q_Btr = radar_comm.Btrdoppler;
% Q_Bmr = radar_comm.Bmrdoppler;
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
for nr = 1:Nr
    for k = 1:K
        q_rnr_k = Qr(:,k,nr);
        Q_rnr_k = diag(q_rnr_k);
        a_k = A(k,:).';
        s_rtr_k = Q_rnr_k*a_k;
        S_rtr(k,:,nr) = s_rtr_k.';
    end
end
S_tr = S_rtr;
R_tr = zeros(K,K,Nr); %Covariance for the radar
for nr = 1:Nr
    eta_rtnr = eta_rtr(nr);
    Sigma_tnr = radar.Sigma(:,:,nr);
    for m = 1:K 
        s_rt_nr_m = S_rtr(m,:,nr).';
        s_tnr_m = S_tr(m,:,nr).';
        for l = 1:K
            s_rt_nr_l = S_rtr(l,:,nr).';
            R_tr(m,l,nr) = eta_rtnr*s_rt_nr_l'*s_rt_nr_m;
            s_tnr_l = S_tr(l,:,nr).';
            R_tr(m,l,nr) = trace(s_tnr_m*s_tnr_l'*Sigma_tnr);
        end
    end
    R_tr(:,:,nr) = R_tr(:,:,nr);
end
cov.S_tr = S_tr;
cov.S_rtr = S_rtr;
cov.target2radar = R_tr;
% radar interference matrices
R_Zr = zeros(K,K,Nr);
R_in_r = zeros(K,K,Nr);
sigma_n = radar.noisepower;
R_r = zeros(K,K,Nr);
for nr = 1:Nr
    R_Zr(:,:,nr) = sigma_n*eye(K,K);
    R_in_r(:,:,nr) = R_Zr(:,:,nr)+ R_Bmr(:,:,nr)+ R_C(:,:,nr);
    R_r(:,:,nr) = R_in_r(:,:,nr)+R_tr(:,:,nr);
end
cov.noise = R_Zr;
cov.inr = R_in_r;
cov.total_r = R_r;
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
        R_rj_k = R_rj(:,:,k);
        R_in_DL{jj,k} = R_MUI_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
        R_dj_k = R_BJ{jj,k};
        R_total_DL{jj,k} = R_dj_k + R_MUI_j_k + R_rj_k + sigma_n*eye(N_j,N_j);
    end
end
cov.in_DL = R_in_DL;
cov.total_DL = R_total_DL;
end

