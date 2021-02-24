function [fdcomm] = Xi_comm_k_UL(fdcomm, k)
Xi_UL_k = 0;
I = fdcomm.UL_num;
for ii = 1:I
    alpha_ui = fdcomm.alpha_UL(ii);
    W_iu_k = fdcomm.UL_weights{ii,k};
    E_iu_k = fdcomm.UL_MMSE_nop{ii,k};
    Xi_UL_k = Xi_UL_k + alpha_ui*real(trace(W_iu_k*E_iu_k));
end
fdcomm.Xi_UL(k) = Xi_UL_k;
end

