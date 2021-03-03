function [radar] = tsp_radar_MMSE_rev(radar, cov)
% radar_MMSE returns the MMSE receiver and Weight Matrice

Nr = radar.Rx;
M = radar.total_Tx;
I_r = zeros(Nr,1); 
K = radar.codelength;
Sigma_t_Nr = radar.Sigma_h_tr;
xi_r = zeros(K,K,Nr);
S_tr = cov.S_tr;
for nr = 1 : Nr
    alpha_nr = radar.alpha_r(nr);
    Sigma_tnr = Sigma_t_Nr(:,:,nr);
    %St_nr = cov.S_tr(:,:,nr);
    R_in_nr = cov.inr(:,:,nr);
    R_t_nr = cov.target2radar(:,:,nr);
    Urnr = Sigma_tnr*S_tr'/(S_tr*Sigma_tnr*S_tr'+R_in_nr);
    R_r_nr = R_in_nr + R_t_nr;
    Ernr_star = Sigma_tnr-Sigma_tnr*S_tr'/(R_r_nr)*S_tr*Sigma_tnr;
    %Wrnr = inv(Ernr_star);
    Wrnr = inv(Ernr_star);
    radar.WMMSE_weights{nr,1} = Wrnr;
    %Ernr = nearestSPD(Sigma_tnr-2*nearestSPD(Urnr*(St_nr*St_nr')*Urnr'))+nearestSPD( Urnr*R_in_nr*Urnr'));
    % Ernr = nearestSPD(Sigma_tnr-2*Urnr*St_nr*Sigma_tnr + nearestSPD(Urnr*(St_nr*St_nr')*Urnr')+nearestSPD(Urnr*R_in_nr*Urnr'));
    % Ernr = Sigma_tnr- Sigma_tnr*St_nr'*Urnr'-Urnr*St_nr*Sigma_tnr+Urnr*St_nr*Sigma_tnr*St_nr'*Urnr'...
%     Ernr = Sigma_tnr- 2*Urnr*S_tr*Sigma_tnr+ Urnr*S_tr*Sigma_tnr*S_tr'*Urnr'...
%         +Urnr*R_in_nr*Urnr';
    radar.WMMSE_RX{nr,1}= Urnr;
    radar.MMSE{nr,1} = Ernr_star;
%     radar.MMSE_nop{nr,1} = Ernr;
    % ee = nearestSPD(eye(M) + nearestSPD(Urnr*R_t_nr*Urnr')/nearestSPD(Urnr*R_in_nr*Urnr'));
    %ee = eye(size(Sigma_t_Nr,2)) + (Urnr*R_t_nr*Urnr')/(Urnr*R_in_nr*Urnr'+1e-8*eye(K*M));
    R_in_tilde = Urnr*R_in_nr*Urnr';
    d1 = eye(size(R_in_tilde), 'logical');
    R_in_tilde(d1) = real(diag(R_in_tilde));
    R_t_tilde = Urnr*R_t_nr*Urnr';
    d2 = eye(size(R_t_tilde), 'logical');
    R_t_tilde(d2) = real(diag(R_in_tilde));
    ee = eye(size(Sigma_t_Nr,2)) + (R_t_tilde)/(R_in_tilde);
    I_r_nr= real(log2(det(ee)));
    I_r(nr)=I_r_nr;
    for k = 1 : K
        urnr_k = Urnr(:,k);
        for m = 1:k
            urnr_m = Urnr(:,m);
            %xi_r(m,k,nr) = alpha_nr*urnr_m.'*Wrnr*conj(urnr_k);
            %xi_r(m,k,nr) =trace(urnr_m.'*Wrnr.'*conj(urnr_k));
            xi_r(m,k,nr) =urnr_k'/(Ernr_star)*urnr_m;
        end
    end
    xi_r(:,:,nr) = alpha_nr*xi_r(:,:,nr);
end
radar.MI_radar = I_r;
radar.xi_r = xi_r;
end

