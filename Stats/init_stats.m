function [stats] = init_stats()
stats=[];
stats.I_total=[];%zeros(Nu,num_drops);
stats.I_total_op = [];
stats.I_UL_op = [];
stats.I_DL_op = [];
stats.I_radar_op = [];
stats.I_max = 0;
stats.I_avg = 0;
 end
