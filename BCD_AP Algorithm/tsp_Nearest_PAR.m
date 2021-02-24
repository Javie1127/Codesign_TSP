function [radar] = tsp_Nearest_PAR(radar)
%function [a_mr_bigstar] = Nearest_PAR(z,rho,c)
%Nearest_PAR returns the nearest a 
%
%% Initialization
Mr = radar.Tx;
A = radar.codematrix;
d = radar.codelength;
for mr = 1:Mr
    rho_mr = radar.gamma_r(mr);
    z_mr = A(:,mr);
    z_nor_mr = normalize(z_mr,'norm',2);
    c = radar.Power(mr);
    sigma_mr = sqrt(c*rho_mr/d);
    M = find(abs(z_nor_mr) == min(abs(z_nor_mr)));
    a_mr_op = zeros(d,1);
    if min(abs(z_nor_mr)) == 0
        for ii = 1:d
            if ismember(ii,M) 
                a_mr_op(ii) = sqrt((c-k*sigma_mr^2)/(d-k));
            else
                a_mr_op(ii) = sigma_mr*exp(1i*angle(z_nor_mr(ii)));
            end
        end
    else
        if numel(unique(z_nor_mr(M))) == 1
            k = d - size(M,1);
        else
            zz = unique(z_nor_mr(M)); % unique elements with the least magnitude
            KK = numel(zz);
            num_KK = zeros(KK,1);
            for kk = 1:KK
                num_KK(kk) = numel(find(z_nor_mr == zz(kk)));
            end
            [KK_min,kk_min] = min(num_KK);
            k = d - (numel(M)- KK_min);
            ind_min = find(z_nor_mr == zz(kk_min));
            M = intersect(M,ind_min);
        end
        gamma=sqrt((c-k*sigma_mr.^2)/norm(z_nor_mr(M))^2);
        if gamma*z_nor_mr(M)>sigma_mr
            for ii = 1:d
                a_mr_op(ii) = sigma_mr*exp(1i*angle(z_nor_mr(ii)));
            end
            %M = find(z_nor_mr(1:d-k) == min(a_mr_star_mag(1:d-k)));
        else
            for ii = 1:d
                if ismember(ii,M) 
                    a_mr_op(ii) = gamma*z_nor_mr(ii);
                else
                    a_mr_op(ii) = sigma_mr*exp(1i*angle(z_nor_mr(ii)));
                end
            end
        end
        radar.codematrix(:,mr) = a_mr_op;
    end
end

