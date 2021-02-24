function [SINR_tot_lin,SINR_tot,Cap_tot_AS, Tput_tot_AS, Rate_oh, SINR_lin_eff, mcs_ues, rank_ues, Rate_rec_oh, detect_outcome]=collect_stats_per_iter(...
    SINR_lin,Cap, Nfreq, Sys_Tput_factor, L, rate_per_layer_cap, RC_output)

SINR_tot_lin=sum(SINR_lin,2)/Nfreq;
SINR_tot=10*log10(SINR_tot_lin);%1-dim vector across all UEs and their SINR averaged over frequency

%Capacity
Cap_tot_AS=sum(Cap,2)/Nfreq;%vector across all UEs and their SINR averaged over frequency
Cap_tot_AS = min(Cap_tot_AS, L*rate_per_layer_cap);

% to generate Tput in unit of Mbps, assume 20MHz(1200RE), 12 data symbols/subframe, pilot overhead 10% 
Tput_tot_AS = Cap_tot_AS*Sys_Tput_factor;

% SINR_lin_eff_dB = 10*log10(SINR_lin_eff);
Rate_oh = RC_output.Rate*0.9; %%% assuming 10% overhead of pilots

Rate_rec_oh = RC_output.Rate_rec*0.9; %%% assuming 10% overhead of pilots

SINR_lin_eff = RC_output.SINR_lin_eff;

mcs_ues = RC_output.mcs_ues;
rank_ues = RC_output.rank_ues;
detect_outcome = RC_output.detect_outcome;

end
