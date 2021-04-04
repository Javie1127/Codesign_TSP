function [radar] = tsp_radar_MMSE_rev(radar, cov)
% radar_MMSE returns the MMSE receiver and Weight Matrice

Nr = radar.Rx;
K = radar.codelength;
Sigma_t_Nr = radar.Sigma_h_tr;
xi_r = zeros(K,K,Nr);
S_tr = cov.S_tr;
R_r = cov.total_r;
%R_in = cov.inr;
for nr = 1 : Nr
    alpha_nr = radar.alpha_r(nr);
    Sigma_tnr = Sigma_t_Nr(:,:,nr);
    %R_in_nr = R_in(:,:,nr);
    R_r_nr = R_r(:,:,nr);
    Urnr = Sigma_tnr*S_tr'/ R_r_nr;
    Ernr_star = nearestSPD(Sigma_tnr-(Sigma_tnr*S_tr'*((R_r_nr\S_tr)*Sigma_tnr)));
%     d_new = eye(size(Ernr_star), 'logical');
%     Ernr_star(d_new) = abs(diag(Ernr_star));
    qq = chol(Ernr_star);
    Wrnr = qq\(qq'\eye(size(Ernr_star,1)));
    d_new = eye(size(Ernr_star), 'logical');
    Wrnr(d_new) = real(diag(Wrnr));
    %Wrnr = nearestSPD(inv_chol(Ernr_star));
    radar.WMMSE_weights{nr,1} = Wrnr;
    radar.WMMSE_RX{nr,1}= Urnr;
    radar.MMSE{nr,1} = Ernr_star;
%     radar.MMSE_nop{nr,1} = Ernr;
    % ee = nearestSPD(eye(M) + nearestSPD(Urnr*R_t_nr*Urnr')/nearestSPD(Urnr*R_in_nr*Urnr'));
    %ee = eye(size(Sigma_t_Nr,2)) + (Urnr*R_t_nr*Urnr')/(Urnr*R_in_nr*Urnr'+1e-8*eye(K*M));
  
    for k = 1 : K
        urnr_k = Urnr(:,k);
        for m = 1:k
            urnr_m = Urnr(:,m);
            %xi_r(m,k,nr) = alpha_nr*urnr_m.'*Wrnr*conj(urnr_k);
            %xi_r(m,k,nr) =trace(urnr_m.'*Wrnr.'*conj(urnr_k));
            xi_r(m,k,nr) =urnr_k'*(Ernr_star\urnr_m);
%             xi_r(m,k,nr) = real(urnr_k'*(Wrnr)*urnr_m);
        end
    end
    xi_r(:,:,nr) = alpha_nr*xi_r(:,:,nr);
end

radar.xi_r = xi_r;
end

