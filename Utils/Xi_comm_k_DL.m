function [DLcomm] = Xi_comm_k_DL(DLcomm, k)
J = DLcomm.DL_num;
Xi_DL_k = 0;
for jj = 1:J
    alpha_dj = DLcomm.alpha_DL(jj);
    W_jd_k = DLcomm.DL_weights{jj,k};
    E_jd_k = DLcomm.DL_MMSE_nop{jj,k};
    Xi_DL_k = Xi_DL_k+alpha_dj*real(trace(W_jd_k*E_jd_k));
end
DLcomm.Xi_DL(k) = Xi_DL_k;
end

