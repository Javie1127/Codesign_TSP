function [inv_A] = inv_SVD(A)
%
[U,S,V] = svd(A);
SS = diag(S);
max_sigma_i = max(SS);
thd = eps*max_sigma_i;
reci_SS = 1./SS;
reci_SS(SS<thd) = 0;
S_inv = diag(reci_SS);
inv_A = V*S_inv*U';
end

