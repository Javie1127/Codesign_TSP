function [fdcomm] = tsp_codesign_BD_DL(fdcomm)
%Block Diagonilzation Precoding
%  
%% Precoding matrices initialization
K = fdcomm.symbol_num_per_frame;
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
    P_DL_BD = cell(J,K);
    P_DL = fdcomm.BS_power;
    d_DL = fdcomm.DLstream_num;
    Mc = fdcomm.BSTx;
    for jj = 1:J
        d_dj = d_DL(jj);
        H_Bj = H_DL{jj};
        H_MUI_j = [];
        for g = 1:J
            H_MUI_j = [H_DL{g}.' H_MUI_j];
        end
        H_MUI_j = H_MUI_j.';
        [~,~,V] = svd(H_MUI_j);
        V_null = V(:,Mc-d_dj+1:end);
        H_ef_j = H_Bj*V_null;
        [~,~,V_ef] = svd(H_ef_j);
        P_dj_ini = V_null*V_ef(:,1:d_dj);
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
            P_DL_BD{jj,k} = P_Bj_ini;
        end
    end
    fdcomm.DLprecoders = P_DL_BD;
end
end

