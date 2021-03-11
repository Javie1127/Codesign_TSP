function [cov] = update_cov_dj_k(fdcomm,radar,cov,radar_comm,P_dj_k_old,jj,k)

%% Update UL related covs
Mc = fdcomm.BSTx;% Number of BS TX antennas
J = fdcomm.DL_num; % Number of DL UEs
H_DL = fdcomm.DLchannels;
D_DL = fdcomm.DLsymbols;
P_dJ = fdcomm.DLprecoders; % cell(J,K*N)
n_Bm = radar_comm.n_Bm;
P_dj_k = P_dJ{jj,k};
R_BJ = cov.DL;

R_MUI_DL = cov.MUI_DL;
R_in_DL  = cov.in_DL;
R_total_DL = cov.total_DL;

for g = 1:J
    HBg = H_DL{g}; %load the UL channel matrix
    if g~= jj
        R_j_MUI_old = R_MUI_DL{g,k};
        R_gj_MUI_old = nearestSPD(HBg*(P_dj_k_old*P_dj_k_old')*HBg');
%         d_old =  eye(size(R_gj_MUI_old), 'logical');
%         R_gj_MUI_old(d_old) = abs(diag(R_gj_MUI_old));
        R_gj_MUI_new = nearestSPD(HBg*(P_dj_k*P_dj_k')*HBg');
%         d_new = eye(size(R_gj_MUI_new), 'logical');
%         R_gj_MUI_new(d_new) = abs(diag(R_gj_MUI_new));
        R_j_MUI_new = R_j_MUI_old - R_gj_MUI_old+R_gj_MUI_new;
        R_MUI_DL{g,k} = R_j_MUI_new;
        R_in_DL{g,k} = R_in_DL{g,k} - R_j_MUI_old+R_j_MUI_new;
        R_total_DL{g,k} = R_total_DL{g,k} - R_j_MUI_old+R_j_MUI_new;
    else
        R_dj_k_old = R_BJ{g,k};
        R_dj_k_new = nearestSPD(HBg*(P_dj_k*P_dj_k')*HBg');
%         d_new = eye(size(R_dj_k_new), 'logical');
%         R_dj_k_new(d_new) = abs(diag(R_dj_k_new));
        R_BJ{g,k} = R_dj_k_new;
        R_total_DL{g,k} = R_total_DL{g,k} - R_dj_k_old + R_dj_k_new;
    end
end
cov.MUI_DL = R_MUI_DL;
cov.DL= R_BJ;
cov.in_DL = R_in_DL;
cov.total_DL = R_total_DL;

%% DL - radar
Nr = radar.Rx;
K = radar.codelength;
S_Bm = cov.S_Bm;
for k = 1:K
    s_dj_k_nBm = S_Bm(:,k);
    P_Bj_k  = P_dJ{jj,k};
    s_dj_k_nBm = P_Bj_k*D_DL{jj}(:,n_Bm,k) + s_dj_k_nBm;
    S_Bm(:,k) = s_dj_k_nBm;
end
R_Bmr = zeros(K,K,Nr);
Sigma_Bm_Nr = radar_comm.Sigma_Bm_Nr;
for nr = 1:Nr
    R_Bm_nr = zeros(K,K);
    for m = 1:K
        S_Bm_m = S_Bm(:,m);
        for l = 1:K
            S_Bm_l = S_Bm(:,l);
            Sigma_Bm_nr_m_l = Sigma_Bm_Nr{m,l,nr};
            R_Bm_nr(m,l) = trace(S_Bm_m*S_Bm_l'*Sigma_Bm_nr_m_l);
        end
    end
    R_Bm_nr = nearestSPD(R_Bm_nr);
    R_Bmr(:,:,nr) = R_Bm_nr;
end
cov.Bmr = R_Bmr;
cov.S_Bm = S_Bm;
%% BS-target-radar
if radar_comm.isCollaborate
    S_Bt = cov.S_Bt;
    S_t_cell = cov.S_t_cell;
    for k = 1:K
        s_dj_k_one_old = S_Bt(:,k);
        s_dj_k_one_new =s_dj_k_one_old- P_dj_k_old*D_DL{jj}(:,1,k)...
            +P_dj_k*D_DL{jj}(:,1,k);
        S_Bt(:,k) = s_dj_k_one_new;
        s_t_k = S_t_cell{k,1}.';
        s_t_k(end-Mc+1:end) = s_dj_k_one_new;
        S_t_cell{k,1} =s_t_k.';
    end
    % https://www.mathworks.com/matlabcentral/answers/46316-sparse-block-diagonal-matrix
    S_tr = blkdiag(S_t_cell{:}); 
    cov.S_tr = S_tr;
    R_Btr = zeros(K,K,Nr);
    Sigma_Bt_Nr = radar_comm.Sigma_Bt_Nr;
    R_rt_r = cov.rtr;
    for nr = 1:Nr
        R_Bt_nr = zeros(K,K);
        for m = 1:K
            S_Bt_m = S_Bt(:,m);
            for l = 1:K
                S_Bt_l = S_Bt(:,l);
                Sigma_Bt_nr_m_l = Sigma_Bt_Nr{m,l,nr};
                R_Bt_nr(m,l) = (trace(S_Bt_m*S_Bt_l'*Sigma_Bt_nr_m_l));
            end
        end
        R_Bt_nr = nearestSPD(R_Bt_nr);
        R_Btr(:,:,nr)=R_Bt_nr;
    end
    R_tr = R_Btr + R_rt_r;
    cov.Btr = R_Btr;
    cov.target2radar = R_tr;
    cov.S_t_cell = S_t_cell;
    cov.S_Bt = S_Bt;
end
%% Update HBB
if fdcomm.UL_num > 0
    I = fdcomm.UL_num;
    H_BB = fdcomm.BBchannel;
    R_BB = cov.B2B;
    R_BB_old = R_BB;
    R_BB_new = nearestSPD(H_BB*(S_Bm(:,k)*S_Bm(:,k)')*H_BB');
%     d_new = eye(size(R_BB_new), 'logical');
%     R_BB_new(d_new) = abs(diag(R_BB_new));
    R_BB{k,1} = R_BB_new;
    cov.B2B = R_BB;
    R_in_UL = cov.in_UL;
    R_total_UL = cov.total_UL;
    for ii = 1:I
        R_in_UL{ii,k} = R_in_UL{ii,k}-R_BB_old{k,1} + R_BB{k,1};
        R_total_UL{ii,k} = R_total_UL{ii,k} - R_BB_old{k,1} + R_BB{k,1};
    end
    cov.in_UL = R_in_UL;
    cov.total_UL = R_total_UL;
end
end