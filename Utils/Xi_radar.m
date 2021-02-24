function [radar] = Xi_radar(radar)
Nr = radar.Rx;
Xi_r = zeros(Nr,1);

for nr = 1:Nr
    alpha_nr = radar.alpha_r(nr);
    W_rnr = radar.WMMSE_weights{nr,1};
    E_rnr = radar.MMSE_nop{nr,1};
    Xi_r_nr = alpha_nr*real(trace(W_rnr*E_rnr));
    Xi_r(nr) = Xi_r_nr;
end
radar.Xi_r = Xi_r;
end

