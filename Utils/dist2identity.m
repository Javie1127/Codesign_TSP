function [dist_A] = dist2identity(A)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
diag_A = real(diag(A));
diag_dist_A = abs(diag_A-ones(size(A,1),1));
diff_A_I = A-eye(size(A,1),1);
d_new = eye(size(diff_A_I), 'logical');
diff_A_I(d_new) = diag_dist_A;
dist_A = norm(diff_A_I,'fro');
end

