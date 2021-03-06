function [fdcomm] = tsp_codesign_Uniform_Precoding_UL(fdcomm)
%Uniform Precoding for UL
%  
%% Precoding matrices initialization
K = fdcomm.symbol_num_per_frame;
%% Precoding matrices initialization
if fdcomm.UL_num>0
    Nc = fdcomm.BSRx;
    I = fdcomm.UL_num;
    P_UL = fdcomm.UL_power;
    P_UL_ini = cell(I,K);
    d_UL = fdcomm.ULstream_num;
    N_UL = fdcomm.UL_UE_Ant;
    for ii = 1:I
        P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
        d_ui    = d_UL(ii);
        N_ui    = N_UL(ii);
        S_ui = eye(N_ui);
        V_ui = eye(d_ui);
        U_ui = [sqrt(P_UL_i)*eye(d_ui);zeros(N_ui-d_ui,d_ui)];
        P_ui_ini = S_ui*U_ui*V_ui;
        for k = 1 : K
            P_UL_ini{ii,k} = P_ui_ini;
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




end

