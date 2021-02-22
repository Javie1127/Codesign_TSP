function [fdcomm, radar, cov] = tsp_precoder_ini(radar,fdcomm,radar_comm,cov)
%%%%-----------------------------------

K = radar.codelength;
%% Precoding matrices initialization
if fdcomm.UL_num>0
    I = fdcomm.UL_num;
    P_UL = fdcomm.UL_power;
    P_UL_ini = cell(I,K);
    d_UL = fdcomm.ULstream_num;
    N_UL = fdcomm.UL_UE_Ant;
    if strcmp(fdcomm.Initializations,'Random')
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
    elseif strcmp(fdcomm.Initializations,'Efficient')
        H_UL = fdcomm.ULchannels;
        for ii = 1:I
            P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
            HiB     = H_UL{ii}; %% UL channel matrix with dimension Mc*Ni
            d_ui    = d_UL(ii);
            [~,~,V] = svd(HiB);
            P_iB_ini = V(:,1:d_ui);
            [U_UL_i,s,V_UL_i] = svd(P_iB_ini);
            s(1:d_ui,:) = diag(sqrt(P_UL_i/d_ui)*ones(d_ui,1));
            %S_sort = diag(s_sort_nonzero);
            P_iB_ini = U_UL_i*s*V_UL_i';
            for k = 1 : K
                P_UL_ini{ii,k} = P_iB_ini;
            end
        end
    end
    fdcomm.ULprecoders = P_UL_ini;
end
if fdcomm.DL_num>0
    H_DL = fdcomm.DLchannels;
    J = fdcomm.DL_num;
    P_DL_ini = cell(J,K);
    P_DL = fdcomm.BS_power;
    d_DL = fdcomm.DLstream_num;
    if strcmp(fdcomm.Initializations,'Random')
        % Initialization approach 2 cholesky decompostion of the diagonal covariance matrix
        Mc = fdcomm.BSTx;
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
    elseif strcmp(fdcomm.Initializations,'Efficient')
        for jj = 1:J
            HBj = H_DL{jj}; %% UL channel matrix with dimension Ni*Mc
            %Mc = size(HBj,2);
            d_dj = d_DL(jj);
            [~,~,V] = svd(HBj);
            P_dj_ini = V(:,1:d_dj);
            [U_DL_j,s,V_DL_j] = svd(P_dj_ini);
            s(1:d_dj,:) = diag(sqrt(P_DL/J/d_dj)*ones(d_dj,1));
            %HBj_H = HBj';
            %P_Bj_ini = HBj_H(:,1:d_Bj);
            %[U_DL_j,~,V_DL_j] = svd(P_Bj_ini);
            %s_DL = diag(S_DL);
            %[~,Idx] = sort(s_DL,'descend');
            %U_sort_j = U_DL(:,Idx);
            %V_sort_j = V_DL(:,Idx);
            %s_sort_nonzero_j = sqrt(P_DL/J/d_Bj)*ones(d_Bj,1);
            %S_sort_j = [diag(s_sort_nonzero_j);zeros(Mc-d_Bj,d_Bj)];
            %     P_Bj_ini = U_DL_j*S_sort_j*V_DL_j';
            P_Bj_ini = U_DL_j*s*V_DL_j';
            for k = 1 : K
                P_DL_ini{jj,k} = P_Bj_ini;
            end
        end
    end
    fdcomm.DLprecoders = P_DL_ini;
end

Mr = radar.Tx;
Pr = radar.Power;
A_ini = zeros(K,Mr);
for mr = 1:Mr
    P_k = Pr(mr);
    A_ini(:,mr) = sqrt(P_k/K);
end
radar.codematrix = A_ini;
S_tr_cell = cell(K,1);
for k = 1:K
    S_tr_cell{k,1} = A_ini(k,:);
end
S_tr = blkdiag(S_tr_cell{:});
cov.S_tr = S_tr;
%% Initializing the covariance matrices
cov = tsp_covmat_rev(fdcomm,radar,radar_comm,cov);
%% Initializing the MMSE matrices
fdcomm = tsp_Comm_MMSE_rev(fdcomm,radar,cov);
radar = tsp_radar_MMSE(radar,cov);
%% Initializing the performance measures
radar = Xi_radar(radar);
for k = 1:K
    fdcomm = Xi_comm_k(fdcomm,k);
end


end
