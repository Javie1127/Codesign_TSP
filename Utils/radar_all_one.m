function [radar] = radar_all_one(radar)

Mr = radar.Tx;
Pr = radar.Power;
K = radar.codelength;
A = zeros(K,Mr);
for mr = 1:Mr
    Pr_mr = Pr(mr);
    A(:,mr) = sqrt(Pr_mr/K)*ones(K,1);
end
radar.codematrix = A;
end

