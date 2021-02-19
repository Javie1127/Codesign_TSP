function [Pd,Pfa] = tsp_radar_detection(fdcomm,radar,radar_comm, M, vu,cov)
%%--------------------%%
%% Test statics deterministic part 
Nr = radar.RX;
Mr = radar.TX;
K = radar.codelength;
Mc = fdcomm.BSTX;
n = radar.CUT_Idx; % CUT index
eta_rtr = radar.channelgain;
eta_Btr = radar_comm.Btrchannelgains;
eta_Bmr = radar_comm.Bmrchannelgains;
if nargin == 5
    cov = tsp_covmat(fdcomm,radar,radar_comm);
end
eta_UL_r = radar_comm.UL2rchannelgains;
N_UL = fdcomm.UL_UE_Ant;
D_UL = fdcomm.ULsymbols;
D_DL = fdcomm.DLsymbols;
sigma0 = radar.noisepower;
R_in_r = cov.inr;
R_tr = cov.target2radar;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
%% Data generation
T_h1 = zeros(M,1); % test statistic for H1
T_h0 = zeros(M,1); % test statistic for H0
for m = 1:M
    T_y_h1 = 0;
    T_y_h0 = 0;
    for nr = 1:Nr
        h_rtr = zeros(Mr,1);
        for mr = 1 : Mr
            h_rtr(mr) = sqrt(eta_rtr(mr,nr)/2)*(randn(1,1)+1i*randn(1,1));
        end
        h_Btr = sqrt(eta_Btr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
        S_tnr = cov.S_tr(:,:,nr);
        h_tr = [h_rtr;h_Btr];
        y_tnr = S_tnr*h_tr;
        sigma_cnr = radar.cluttercov(:,:,nr);
        y_cnr_re = mvnrnd(zeros(K,1),sigma_cnr/2).';
        y_cnr_im = mvnrnd(zeros(K,1),sigma_cnr/2).';
        y_cnr = y_cnr_re + 1i*y_cnr_im;
        y_Unr = zeros(K,1);
        for k = 1:K
            y_Unr_k = 0;
            for ii = 1:I
                H_i_nr = sqrt(eta_UL_r(ii,nr)/2)*(randn(N_UL(ii),1)+1i*randn(N_UL(ii),1));
                P_ui_k = fdcomm.ULprecoders{ii,k};
                d_ui_k_n = D_UL{ii}(:,1+n,k);
                y_Unr_k = H_i_nr'*(P_ui_k * d_ui_k_n) + y_Unr_k;
            end
            y_Unr(k) = y_Unr_k;
        end
        y_Bmnr = zeros(K,1);
        for k = 1:K
            y_Bmnr_k = 0;
            for jj = 1:J
                H_Bm_nr = sqrt(eta_Bmr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
                P_dj_k = fdcomm.DLprecoders{jj,k};
                d_dj_k = D_DL{jj}(:,n+1,k);
                y_Bmnr_k = H_Bm_nr' * (P_dj_k*d_dj_k) + y_Bmnr_k;
            end
            y_Bmnr(k) = y_Bmnr_k;
        end
        z_nr = sqrt(sigma0/2).*(randn(K,1)+1i*randn(K,1));
        y_rnr_h1 = y_tnr+y_cnr+y_Bmnr+y_Unr+z_nr;
        y_rnr_h0 = y_cnr+y_Bmnr+y_Unr+z_nr;
        Mnr = y_rnr_h0'*y_rnr_h0/K;
        %Mnr = R_in_r(:,:,nr);
        y_bar_rnr_h1 = Mnr^(-1/2)*y_rnr_h1;
        y_bar_rnr_h0 = Mnr^(-1/2)*y_rnr_h0;
        % R_t_nr = R_tr(:,:,nr) ;
        R_t_nr = y_tnr'*y_tnr/K;
        G_nr = Mnr^(-1/2)*R_t_nr*Mnr^(-1/2);
        [V_nr,Lambda_nr] = eig(G_nr);
        T_y_nr_h1 = y_bar_rnr_h1'*V_nr/(Lambda_nr^(-1)+eye(K))*V_nr'*y_bar_rnr_h1;
        T_y_nr_h0 = y_bar_rnr_h0'*V_nr/(Lambda_nr^(-1)+eye(K))*V_nr'*y_bar_rnr_h0;
        T_y_h1 = T_y_h1 + T_y_nr_h1;
        T_y_h0 = T_y_h0 + T_y_nr_h0;
    end
    T_h1(m) = T_y_h1;
    T_h0(m) = T_y_h0;
end
num_h1 = find(T_h1>vu); % detection
num_h0 = find(T_h0>vu);% false alarm
Pd = length(num_h1)/M;
Pfa = length(num_h0)/M;

end
