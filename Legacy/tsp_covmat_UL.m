function [cov] = tsp_covmat_UL(fdcomm,radar,radar_comm)

%% radar only covariances
Mr = radar.TX; % Number of radar TX antennas
Nr = radar.RX; % Number of radar RX antennas
K = radar.codelength; % The length of the radar code; or the number of PRIs
n = radar.CUT_Idx; % CUT index
%H_rtr                       = radar.channel;
Qr = radar.doppler; %Mr*K*Nr radar temporal steering vector
R_C = radar.cluttercov;
A = radar.codematrix;
sigma_0 = radar.noisepower;
%% UL comm related covariance matrices
Nc = fdcomm.BSRX;% Number of BS RX antennas
I = fdcomm.UL_num; % Number of UL UEs
P_UI = fdcomm.ULprecoders; % cell(I,k*N)
D_UL = fdcomm.ULsymbols;
H_UL = fdcomm.ULchannels;
% eta_BB = fdcomm.BBchannelgain;
R_IB = cell(I,K);   % All the UL covariance matrices
R_MUI_UL = cell(I,K);
R_sum_UL = cell(K,1);
%% radar-comm 
eta_UL_r = radar_comm.UL2rchannelgains;
Q_Ir = radar_comm.ULrdoppler;
%% Covariance matrix at the nr radar RX 
% covariance matrix R_t_r 
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
    Sigma_tnr = radar.Sigma(:,:,nr);
    for m = 1:K 
        s_tnr_m = S_tr(m,:,nr).';
        for l = 1:K
            s_tnr_l = S_tr(l,:,nr).';
            R_tr(m,l,nr) = trace(s_tnr_m*s_tnr_l'*Sigma_tnr);
        end
    end
    R_tr(:,:,nr) = nearestSPD(R_tr(:,:,nr));
end
cov.S_tr = S_tr;
cov.S_rtr = S_rtr;
cov.target2radar = R_tr;
% UL UEs to radar RX
%H_UL_r = radar_comm.UL2radarchannels;
R_Ur = zeros(K,K,Nr);
for nr = 1:Nr
    for m= 1:K
        for l = 1:K
            R_U_nr_temp = 0;
            for ii = 1:I
                P_ui_m = P_UI{ii,m};
                d_ui_m = D_UL{ii}(:,1+n,m);
                s_inr_m = P_ui_m*d_ui_m;
                f_i_nr = log(Q_Ir(2,nr,ii))/(2*pi*1i);
                P_ui_l = P_UI{ii,l};
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
    R_in_r(:,:,nr) = R_Zr(:,:,nr)+ R_Ur(:,:,nr)+ R_C(:,:,nr);
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
        PiB_k = P_UI{ii,k};
        R_IB{ii,k} = H_iB*(PiB_k*PiB_k')*H_iB';
        R_sum = R_IB{ii,k}+R_sum;
    end
    R_sum = R_sum;
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
        R_rB_k = R_rB(:,:,k);
        R_iB_k = R_IB{ii,k};
        R_in_UL{ii,k} = R_MUI_i_k + R_rB_k + R_ZB;
        R_total_UL{ii,k} = R_in_UL{ii,k} + R_iB_k;
    end
end
cov.in_UL = R_in_UL;
cov.total_UL = R_total_UL;
end

