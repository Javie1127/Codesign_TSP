function [radar] = randomcode(radar)
Mr = radar.Tx;
Pr = radar.Power;
K = radar.codelength;
A = zeros(K,Mr);
A_ini = randn(K) + 1i*randn(K);
[U_ini,~] = qr(A_ini);
for mr = 1:Mr
    Pr_mr = Pr(mr);
    ii = mod(mr,K);
    if ii == 0
        ii = K;
    end
    a_mr = U_ini(:,ii);
    A(:,mr) = sqrt(Pr_mr)*a_mr;
end
radar.codematrix = A;
end

