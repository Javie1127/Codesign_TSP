function [A] = randomcode(Pr,Mr,K)
A = zeros(K,Mr);
A_ini = randn(K) + 1i*randn(K);
[U_ini,~] = qr(A_ini);
for mr = 1:Mr
    Pr_mr = Pr(mr);
    a_mr = U_ini(:,mr);
    A(:,mr) = sqrt(Pr_mr/K)*a_mr;
end
end

