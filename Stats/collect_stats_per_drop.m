function stats=collect_stats_per_drop(stats, SNR_asoc,SNR,PLs_asoc,...
    Prx_asoc,Cap_tot,num_iter,SINR_tot_lin,~, Tput_tot, Rate_oh_tot, SINR_lin_eff_tot, Rate_rec_oh_tot, updates_overiters)

    stat_cap=(1/num_iter)*sum(Cap_tot,2); %vector of UE capacity avaeraged over all iterations for a given drop            
    stat_sinr=10*log10((1/num_iter)*sum(SINR_tot_lin,2));
    
    % averaged Tput over iterations (averaged over channel fades)
    stat_tput = sum(sum(Tput_tot))/num_iter;%changed to sum over iterations
    
    stats.snr_asoc=[stats.snr_asoc;SNR_asoc];
    stats.snr=[stats.snr;SNR];
    stats.pls=[stats.pls;PLs_asoc];
    stats.prx_asoc=[stats.prx_asoc;Prx_asoc];
    stats.cap=[stats.cap;stat_cap];
    stats.sinr=[stats.sinr;stat_sinr];
    stats.tput=[stats.tput;stat_tput];
    
    rate_ave = (1/num_iter)*sum(Rate_oh_tot,2);
    stats.rate = [stats.rate; rate_ave.'];
    
    rate_rec_ave = (1/num_iter)*sum(Rate_rec_oh_tot,2);
    stats.rate_rec = [stats.rate_rec; rate_rec_ave.'];

    sinr_ave = sum(SINR_lin_eff_tot,2)./sum(SINR_lin_eff_tot>0,2);
    stats.sinrdB_ave = [stats.sinrdB_ave; 10*log10(sinr_ave).'];
    
    stats.sinrdB_periter = [stats.sinrdB_periter; 10*log10(SINR_lin_eff_tot)];
    
    if isfield(updates_overiters,'traffic_model')
        % bursty traffic model related stats
        traffic_model = updates_overiters.traffic_model;
        Nu = length(traffic_model.UE_buffer);
        UPT_drop = zeros(Nu,1);
        latency_drop = zeros(Nu,1);
        success_rate = zeros(Nu,1); % collect successful transmission rate 
        for nu = 1:Nu
            % calculate the packet throughput (TR 36.889 A 1.1 performance
            % metrics) for each UE
            pkt_size_series_nu = traffic_model.pkt_size{nu};
            pkt_size_track_nu = traffic_model.pkt_size_track{nu};
            pkt_nack_nu = traffic_model.ARQ_buffer{nu};
            pkt_stamps_flag_nu = traffic_model.pkt_stamp_flags{nu};
            pkts_transmitted_nu = pkt_size_series_nu-pkt_size_track_nu-pkt_nack_nu;
            pkt_stamps_nu = traffic_model.pkt_stamps{nu};
            pkt_tput_nu = pkts_transmitted_nu./(pkt_stamps_nu*1e-3)/8; % calculate the throughput in bytes/s
            %pkt_tput_nu(isnan(pkt_stamps_flag_nu))=0; % dropped file to
            %have zero throughput 
            pkt_stamps_nu_nondrop = pkt_stamps_nu(~isnan(pkt_stamps_flag_nu));
            pkt_stamps_nu_success = pkt_stamps_nu_nondrop(pkt_stamps_nu_nondrop);
            UPT_drop(nu) = mean(pkt_tput_nu);
            latency_drop(nu) = mean(pkt_stamps_nu_success); % calculated latency out of the successfully transmitted packets
            success_rate(nu) = length(pkt_stamps_nu_success)/length(pkt_stamps_nu);
        end
        stats.UPT = [stats.UPT UPT_drop];
        stats.latency = [stats.latency latency_drop];
        stats.lambda = [stats.lambda traffic_model.arrival_rate.'];
        stats.S = [stats.S traffic_model.file_size.'];
        stats.success_rate = [stats.success_rate success_rate.'];
    end
    
end

