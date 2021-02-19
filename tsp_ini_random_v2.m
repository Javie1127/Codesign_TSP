function [fdcomm, radar, cov] = tsp_ini_random_v2(radar,fdcomm,radar_comm)
%%%%-----------------------------------
%------- SNR.SNR_rtr: Mr*Nr real matrix
%------- SNR.SNR_DL:
%% Precoding matrices initialization
I = fdcomm.UL_num;
K = radar.codelength;
P_UL = fdcomm.ULpower;
P_UL_ini = cell(I,K); 
% %% Initialization approach 1 right sigular matrix 
% for ii = 1:I
%     P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
%     HiB     = H_UL{ii}; %% UL channel matrix with dimension Mc*Ni
%     d_ui    = d_UL(ii);
%     [~,~,V] = svd(HiB);
%     P_iB_ini = V(:,1:d_ui);
%     [U_UL_i,s,V_UL_i] = svd(P_iB_ini);
%     s(1:d_ui,:) = diag(sqrt(P_UL_i/d_ui)*ones(d_ui,1));
%     %S_sort = diag(s_sort_nonzero);
%     P_iB_ini = U_UL_i*s*V_UL_i';
%     for k = 1 : K
%         P_UL_ini{ii,k} = P_iB_ini;
%     end
% end
J = fdcomm.DL_num;
P_DL_ini = cell(J,K);
P_DL = fdcomm.DLpower;
% 
% for jj = 1:J
%     HBj = H_DL{jj}; %% UL channel matrix with dimension Ni*Mc
%     %Mc = size(HBj,2);
%     d_dj = d_DL(jj);
%     [~,~,V] = svd(HBj);
%     P_dj_ini = V(:,1:d_dj);
%     [U_DL_j,s,V_DL_j] = svd(P_dj_ini);
%     s(1:d_dj,:) = diag(sqrt(P_DL/J/d_dj)*ones(d_dj,1));
%     %HBj_H = HBj';
%     %P_Bj_ini = HBj_H(:,1:d_Bj);
%     %[U_DL_j,~,V_DL_j] = svd(P_Bj_ini);
%     %s_DL = diag(S_DL);
%     %[~,Idx] = sort(s_DL,'descend');
%     %U_sort_j = U_DL(:,Idx);
%     %V_sort_j = V_DL(:,Idx);
%     %s_sort_nonzero_j = sqrt(P_DL/J/d_Bj)*ones(d_Bj,1);
%     %S_sort_j = [diag(s_sort_nonzero_j);zeros(Mc-d_Bj,d_Bj)];
% %     P_Bj_ini = U_DL_j*S_sort_j*V_DL_j';
%     P_Bj_ini = U_DL_j*s*V_DL_j';
%     for k = 1 : K
%         P_DL_ini{jj,k} = P_Bj_ini;
%     end
% end
%% Initialization approach 2 cholesky decompostion of the diagonal covariance matrix
d_UL = fdcomm.ULstream_num;
N_UL = fdcomm.UL_UE_Ant;
d_DL = fdcomm.DLstream_num;
for ii = 1:I
    P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
    d_ui    = d_UL(ii);
    N_ui    = N_UL(ii);
    Sigma_ui = [sqrt(P_UL_i/N_ui)*eye(d_ui);zeros(N_ui-d_ui,d_ui)];
    A1 = randn(N_ui) + 1i*randn(N_ui);
    [S_ui,~] = qr(A1);
    A2 = randn(d_ui) + 1i*randn(d_ui);
    [U_ui,~] = qr(A2);
    P_ui_ini = S_ui*Sigma_ui*U_ui';
    for k = 1 : K
        P_UL_ini{ii,k} = P_ui_ini;
    end
end
Mc = fdcomm.BSTX;
for jj = 1:J
    d_dj = d_DL(jj);
    Sigma_dj = [sqrt(P_DL/J)*eye(d_dj);zeros(Mc-d_dj,d_dj)];
    A1 = randn(Mc) + 1i*randn(Mc);
    [S_dj,~] = qr(A1);
    A2 = randn(d_dj) + 1i*randn(d_dj);
    [U_dj,~] = qr(A2);
    P_dj_ini = S_dj*Sigma_dj*U_dj';
    for k = 1 : K
        P_DL_ini{jj,k} = P_dj_ini;
    end
end
Mr = radar.TX;
Pr = radar.Pr;
A_ini = zeros(K,Mr);
for mr = 1:Mr
    P_k = Pr(mr);
    A_ini(:,mr) = sqrt(P_k/K);
end
radar.codematrix = A_ini;
fdcomm.ULprecoders = P_UL_ini;
fdcomm.DLprecoders = P_DL_ini;
%% Initializing the covariance matrices
cov = tsp_covmat(fdcomm,radar,radar_comm);
% %% Initializing the MMSE matrices
% fdcomm = tsp_Comm_MMSE(fdcomm,radar,cov);
% radar = tsp_radar_MMSE(radar,cov);
% %% Initializing the performance measures
% radar = Xi_radar(radar);
% for k = 1:K
%     fdcomm = Xi_comm_k(fdcomm,k);
% end


end
