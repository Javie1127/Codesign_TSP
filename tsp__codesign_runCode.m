function runCode(configureXML)
% if ~isdeployed
%     if isunix
%         addpath(genpath([pwd,'/Utils']));
%         updateMatlabPath()
%     else
%         %have to run from Matlab folder
%         addpath(genpath([pwd, '/Utils']));
%         updateMatlabPath()
%     end
if ~isdeployed
    if isunix
        addpath(genpath([pwd,'/Utils']));
        updateMatlabPath()
    else
        %have to run from Matlab folder
        addpath(genpath([pwd, '/Utils']));
        updateMatlabPath()
    end
end


%end

%read XML file
if nargin==0
    clear;
%     configureXML = 'sanity_checks/case_ul_003/system_cfg.xml';
    configureXML = 'sanity_checks/case_021/system_cfg.xml';
%     configureXML = 'system_cfg';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
systemCfg = set_default_config_parameters(systemCfg);
XCoMP_config_readXML;

% open files for reading writing channels
fids_struct = open_files(storeChan);

% initialize some variables before the loop over drops
[stats, RC_output] = init_stats();
percOffload=zeros(num_drops,1);

if SC_outdoor.isTrackMacroInterf
    stats_hetro=Hetro_INR_init();
end




%Initializaions
if systemCfg.gNB.Rate_Control.RankSelEnable == 1
    sinr_offset_ues_iters = ones(Nu, systemCfg.gNB.Rate_Control.MaximumRank, num_iter);
else
    sinr_offset_ues_iters = ones(Nu, 1, num_iter);
end
sinr_offset_ues_drops = cell(0);
detection_ues_iters = ones(Nu, num_iter, num_drops);
% SINR_delta_test = ones(Nu, num_iter, num_drops);
mcs_ues_tot = nan(Nu, num_iter);
rank_ues_tot = zeros(Nu, num_iter);
detect_outcome_tot = nan(Nu, num_iter);
Rate_oh_tot = zeros(Nu, num_iter);
SINR_lin_eff_tot = zeros(Nu, num_iter);
Rate_rec_oh_tot = zeros(Nu, num_iter);


%%% The loop over drops
for d=1:num_drops
    Nu=dims.Nu;
    Nb=dims.Nb;
    %get a layout, path loss, and gNB to UE asoc
    if systemCfg.ChannelSettings.Quadriga_used
        [coordg_x,coordg_y,coordu_x,coordu_y,coordu_z, TimeStamp, ang_z, ang_y,ang_x]=...
            createLO_qd(LO,layout_len, layout_wid, Lo,Wo,systemCfg.UE.UEHeight, num_gNB_rows, Nb,Nu,Nx,...
            isLOPlot,fig_idx,offset_x, minDistance, deply_SC, groupRaduis, num_group, TrackCoordinatesFromFile, secNumber);%put isPlot=0
        
        % number of gNB antenna in a panel is NtxNt2
        [PLs, LOS, gnbUEBeam_asoc_org, ch_coeff, BlockedAPs, gnbUEBeam_asoc_org_diffQAll] = getChannelMatrix_qd_All(gnb_arry,Nr,Nt,Nt2,Nu,Nb,BW,...
            coordg_x,coordg_y,systemCfg.NetworkLayout.Height,coordu_x,coordu_y, coordu_z,ang_z, ang_y,ang_x, TimeStamp ,fade_gran, fc, Fixed_codebook,...
            layout_len, layout_wid,UE_Pattern,UE_orientation, QuarterBlokageEnabled,chanPerturb, ...
            TrackCoordinatesFromFile, TopBeamsDiffQuadrants);  %PLs_pat is the PLs plus power of antenna pattern
        
        [PLs_asoc,idx_asoc_ue,idx_asoc_gnb,Prx_asoc,Prx, Prx_UL]=ue2gnb_asoc_qd_All(PLs,Nu,Nb,PtxG, gnbUEBeam_asoc_org, ...
            TrackCoordinatesFromFile,isUE_Assoc_Enabled); % 
    else
        [coordg_x,coordg_y,coordu_x,coordu_y,coordx_x,coordx_y,Nb_gNB]=...
            createLO(LO,layout_len, layout_wid, Lo,Wo, num_gNB_rows, Nb,Nu,Nx,isLOPlot,fig_idx,offset_x, ...
            systemCfg.NetworkLayout);
        
        [PL,PLs,LOS]=getPLMatrix(deployment_model,fc,coordg_x,coordg_y,coordu_x,...
                coordu_y,height_struct.gnb,height_struct.ue,...
                bypass_shadowing, MCL, systemCfg.NetworkLayout.force_NLOS, ...
                systemCfg.NetworkLayout, systemCfg.gNB.AntennaPattern.enabled, systemCfg.gNB.AntennaPattern, systemCfg.NetworkLayout.Wraparound,Nu, 'macro_cell');
        [PLs_asoc,idx_asoc_ue,idx_asoc_gnb,Prx_asoc,Prx,Prx_UL]=ue2gnb_asoc(PLs,Nu,Nb_gNB,PtxG, PtxUEG,isUE_Assoc_Enabled,systemCfg.UL);
    end
    
    %%% if off-loading
    if Nx~=0
        if rrue_type=="snr_shp"
            [PL,PLs,LOS,PLs_asoc,idx_asoc_ue,idx_asoc_gnb,...
                Prx_asoc,Prx,percOffload(d),Nu,~]=RRUE_ss(deployment_model,...
                fc,coordg_x,coordg_y,coordx_x,...
                coordx_y,coordu_x,coordu_y, height,bypass_shadowing,MCL,Nx,Nb_gNB,Nu,Nofl,PtxG,PrxUEG,PtxGx,...
                PLs_asoc,PL,PLs,LOS);
        elseif rrue_type=="dim_shp"
            [PL,PLs,LOS,PLs_asoc,idx_asoc_ue,idx_asoc_gnb,Prx_asoc,...
                Prx,Nball]= RRUE_dds(deployment_model,fc,coordx_x,...
                coordx_y,coordu_x,coordu_y, height,bypass_shadowing,MCL,Nb_gNB,Nu,PtxG,PtxUEG,...
                PL,PLs,LOS,systemCfg.NetworkLayout.force_NLOS);
        elseif rrue_type=="dim_shp_rept"
            [PLxg,PLsxg,LOSxg,PLxu,PLsxu,LOSxu,PLsxg_asoc,idx_asoc_x2g,...
                idx_asoc_gnbx,Prx_asoc_g2x,Prx_g2x,PLsxu_asoc,...
                idx_asoc_u2x,idx_asoc_xasgnb,Prx_asoc_x2u,Prx_x2u,...
                SNR_asoc_g2x,SNR_g2x,SNR_lin_g2x,...
                SNR_asoc_x2u,SNR_x2u,SNR_lin_x2u,SNR_fac_x]= RRUE_ds_rpt(deployment_model,...
                fc,coordg_x,coordg_y,coordx_x,...
                coordx_y,coordu_x,coordu_y, height,bypass_shadowing,...
                MCL,Nx,Nb_gNB,Nu,PtxG,PtxGx,noise_floor,SNR_cap);
            %%%%%%%%%%%%%%%%%%ADDED for SC
        elseif rrue_type=="outdoor_sc"
            forceNLOS_sc=0;
            [PL,PLs,LOS,PLs_asoc,idx_asoc_ue,idx_asoc_gnb,Prx_asoc,...
                Prx,Nb, hetro_intf]= RRUE_sc(systemCfg.NetworkLayout.RelayNodes.LayoutModel_x,fc,coordx_x,...
                coordx_y,coordu_x,coordu_y, height_struct,bypass_shadowing,MCL_struct.sc,Nu,PtxG,PtxUEG,PtxGx,...
                PL,PLs,LOS,forceNLOS_sc, systemCfg.NetworkLayout, systemCfg.NetworkLayout.Wraparound, systemCfg.gNB,systemCfg.HetNet.Num_sector);
            if  SC_outdoor.isTrackMacroInterf
                [stats_hetro]=Hetro_update_interf(hetro_intf,stats_hetro,noise_floor);
            end
        else
            assert('Nx~=0 and RRUE type not defined');
        end
    end
    
    %%% init variables per drop
    ffinfo_cluster = initFastFadeInfoForClustering(systemCfg, PLs_asoc, Nr, Nr_ob, L, Nt, Nfreq);
    Cap_tot=zeros(Nu,num_iter);
    SINR_tot=zeros(Nu,num_iter);
    SINR_tot_lin=zeros(Nu,num_iter);
    Tput_tot=zeros(Nu,num_iter);
    Bits_tot = zeros(Nu,num_iter);
    idx_AS_ue_log=cell(num_iter,1);
    idx_AS_gnb_log=cell(num_iter,1);
    updates_overiters.ue_priority = ones(Nu,1);
    Indx_ConsRetx = zeros(Nu,1);
    prevServingAP = zeros(Nu,1);
    reTX_TabuList = cell(Nu,1); 
    rank_ues = zeros(1, Nu);
    mcs_ues = nan(Nu, 1);
    SINR_lin_eff = zeros(Nu, 1);
    detect_outcome = nan(Nu, 1);
    Rate = zeros(Nu,1);
    Rate_rec = zeros(Nu,1);
%     HARQ_count = zeros(Nu,1);  %%% for test
%     SINRoffset_accu_HARQ = ones(Nu,2);  %%% for test
    if systemCfg.gNB.Rate_Control.RankSelEnable == 1
        updates_overiters.sinr_offset_ues = ones(Nu, systemCfg.gNB.Rate_Control.MaximumRank);  %%% linear SINR
    else
        updates_overiters.sinr_offset_ues = ones(Nu, systemCfg.UE.NumMIMOLayers);  %%% linear SINR
    end
    
    INR_macro=zeros(Nu,1);
    inrmacro=zeros(Nu,1);
    INR_sc2macro=zeros(Nu,1);
    inrsc2macro=zeros(Nu,1);

    %compute SNR
    [SNR_asoc,SNR,SNR_lin,SNR_UL,SNR_UL_lin]=computeSNR(Prx_asoc,Prx,noise_floor,SNR_cap,Prx_UL);
    
    %%% Calculate MIMO antenna correlation matrix
%     [corr_BS_matrix_all, corr_UE_matrix_all] = antCorrMatrix(systemCfg.ChannelSettings, Nr,Nt,Nu,Nb_gNB);
    
    %Clusterign algrithms based on DL received SNR
    if clustering_type == "single"
        num_AS = 1;
        idx_AS_gnb = cell(1, num_AS);
        idx_AS_ue = cell(1, num_AS);
        idx_AS_gnb{1} = 1:Nb;
        idx_AS_ue{1} = 1:Nu;
        num_AS_gnb_used = Nb;
    elseif any(contains({'capBased_AS', 'capBased_AS_init_off'}, clustering_type))
        if systemCfg.Scheduler.Enabled ~= 1
            % if scheduler is enabled clustering changes every iteration
            [idx_AS_ue,idx_AS_gnb,~,num_AS]=capBased_AS(idx_asoc_ue,idx_asoc_gnb,SNR_lin,Nu,Nb,Nr_ob,Nt,...
                overdim_factor,AS_metric,AS_size_max,O_complexity,lkg_exponent,max_lkg_allowed, updates_overiters.ue_priority,...
                systemCfg.gNB.BBU.Cell_Association,L);
        end
    elseif clustering_type == "fixed_geographic_AS"
        if systemCfg.Scheduler.Enabled ~= 1
            [idx_AS_ue, idx_AS_gnb, num_AS] = fixed_geographic_AS(idx_asoc_ue, ...
                idx_asoc_gnb, coordg_x, coordg_y, AS_size_max, SNR_lin, systemCfg.gNB.BBU.Cell_Association, ...
                updates_overiters.ue_priority);
            num_AS_gnb_used = ones(1, num_AS);
        end
    elseif clustering_type=="hetrogeneous"
        [idx_AS_ue, idx_AS_gnb, num_AS, num_AS_gnb_used] = fixed_outdoor_clusters (idx_asoc_ue, ...
            idx_asoc_gnb,Nx, SC_outdoor.hotspot_numHotspots,SC_outdoor.reuse1, 0, 0);
        idx_AS_ue_orig=idx_AS_ue;
        idx_AS_gnb_orig=idx_AS_gnb;
    elseif any(contains({'ITUGA', 'on_off','DAS','extra_user1','extra_user2'}, clustering_type))
    else
        assert(false, join(['unsupported clustering_type: ', clustering_type]))
    end
    
   if systemCfg.ChannelSettings.Quadriga_used
       tm = initTrafficModel(Nu, num_iter, st80211NR,st80211ad, Slicing);
      [Ht_all, Ht_snr_all]=storedChannelMatrix_qd_perPos(ch_coeff, Nr,Nt, Nt2,Nu,Nb,BW,fade_gran, Prx, PtxG,noise_floor);
       Htot_rc= 0; Htot_h_rc= 0;
       Ht_snr_rc=  0; Ht_snr_h_rc= 0;
       sim_time =0;
       updates_overiters.ue_priority =-Inf*ones(Nu,1);%(sim_time-tm.tdeadline(:,1));
       if TrackCoordinatesFromFile
           sim_time2 = TimeStamp(1)*1e-7;  % to track time stamp
       else
           sim_time2 = 0;
       end
       p = 1; % index for channel update, UE position in the track
       pa=1;   %indx for Beam aging
       dn = zeros(Nu,1); %backoff for CQI
       p_store =[];
       pa_store =[];
   else
        ue_tput_mbps= mean(Tput_tot,2);
        updates_overiters.ue_priority = 1./(1e-6*rand(size(ue_tput_mbps))+ue_tput_mbps);
   end
   
   %%% Add traffic model
   if systemCfg.Scheduler.Enabled && systemCfg.Traffic_model.Enabled
       %% Generate packet arrival times for all the UEs based on a certain FTP model
       %assert(systemCfg.Scheduler.Enabled==0,'Scheduler has to be enabled to deploy Traffic model')
       updates_overiters = traffic_generation(systemCfg,idx_asoc_ue,updates_overiters);
   end
   
   %%% The loop over time slots   
   for iter = 1:num_iter 
       displayProgress(d, num_drops, iter, num_iter, configureXML, isdeployed);
       if systemCfg.ChannelSettings.Quadriga_used
           [p,pa,gnbUEBeam_asoc, gnbUEBeam_asoc_org_diffQ] =BeamUpdateandAging(gnbUEBeam_asoc_org,beamAging, BAfactor, iter, ...
               p,pa, TrackCoordinatesFromFile, TimeStamp,sim_time2 , gnbUEBeam_asoc_org_diffQAll);
           p_store =[p_store  p];
           pa_store =[pa_store pa];
       end
 
       if systemCfg.Scheduler.Enabled
           if ~systemCfg.ChannelSettings.Quadriga_used
               % update UE priority
               ue_tput_mbps= mean(Tput_tot,2);
               updates_overiters.ue_priority = 1./(1e-6*rand(size(ue_tput_mbps))+ue_tput_mbps);
               % EGoS fairness for now, dither the priority to break ties, as
               % ties may impact pruning in some clustering algorithms
           end
           assert(all(isreal(updates_overiters.ue_priority)));
           assert(~any(isnan(updates_overiters.ue_priority)));
       end
      
        %get fading channel samples, weighted and not with SNR
        if ~storeChan.isReadChannel
            if systemCfg.ChannelSettings.Quadriga_used
                [Htot,Htot_h,Ht_snr,Ht_snr_h]=getChannelMatrix_qd_perIter(ch_coeff(:,:,:,p),p,Nr,Nt, Nt2,Nu,Nb,BW,...
                    fade_gran, fc, SNR_lin(:,:,p), Prx(:,:,p), PtxG,noise_floor,chanPerturb );

            else
                if isUL
                    snrlin=SNR_UL_lin;
                    
                else
                    snrlin=SNR_lin;
                end
                if isUL && (UL_Clustering== "ON")
                  snrlin_cl=SNR_UL_lin; 
                else
                  snrlin_cl=SNR_lin;  
                end
                [Htot,Htot_h,Ht_snr,Ht_snr_h]=getChannelMatrix(fadingModel,Nr,Nt,Nu,Nb_gNB,BW,...
                    fade_gran,LOS(:,1:Nb_gNB),Krician,chanPerturb,phase_err,mag_err,...
                    snrlin(:,1:Nb_gNB),isChan_err, [], [], [], [], systemCfg.ChannelSettings); %cell of a size equal number of fades simulated. Each cell has Nu x Nb cells each has Nr x Nt channel matrix
           end
           if and(Nx>0,strcmp(rrue_type,"dim_shp_rept"))
               [Htot,Htot_h,Ht_snr,Ht_snr_h]=getChannelMatrix_dsrpt(Htot,...
                    Htot_h,Ht_snr,Ht_snr_h,...
                    LOSxg,LOSxu,SNR_lin_g2x,SNR_lin_x2u,SNR_fac_x,...
                    Nr,Nt,Nu,Nb_gNB,Nx,Nax,Nfreq,BW,...
                    fade_gran,Krician,chanPerturb,phase_err,mag_err,isChan_err);

            elseif and(Nx * sum(systemCfg.NetworkLayout.SC_outdoor.hotspot_numHotspots)>0, strcmp(rrue_type,"outdoor_sc"))   %%% get the channels for outdoor small cells and concatenate them with gNB channels
                LOS_scs = LOS(:,Nb_gNB+1:Nb);
                snrlin_scs = snrlin(:,Nb_gNB+1:Nb);
                [Htot,Htot_h,Ht_snr,Ht_snr_h]=getChannelMatrix(fadingModel,Nr,Nax,Nu,(Nb-Nb_gNB),BW,...
                    fade_gran,LOS_scs,Krician,chanPerturb,phase_err,mag_err,...
                    snrlin_scs,isChan_err,Htot,Htot_h,Ht_snr,Ht_snr_h,[]);
                
            end
            if storeChan.isSaveChannel
                saveChannelsToFile(fids_struct, Htot, Htot_h, Ht_snr, Ht_snr_h, dims);
            end
        else
            [Htot, Htot_h, Ht_snr, Ht_snr_h]= readChannelsFromFile(fids_struct, dims);
        end
        
      if systemCfg.ChannelSettings.Quadriga_used
         % Consider only UEs with data on the Queue also update Priority accordingly
        [tm, updates_overiters.ue_priority, UEsWithData, idx_asoc_ue_iter,idx_asoc_gnb]= ...
            SchdUEswithDataontheQueue_v2(gnbUEBeam_asoc,sim_time,TimeSlot, updates_overiters.ue_priority, iter, Nu, Nb, tm);
       end
       
       if systemCfg.gNB.BBU.Cell_Association.use_fast_fade_info_for_clustering
            ffinfo_cluster.H_snr=Ht_snr_h;
            ffinfo_cluster.H=Htot_h;
       end
       %run one iteration of code get SINR, cap
       num_layers_for_wsr_approx = systemCfg.UE.NumMIMOLayers;
       if systemCfg.Scheduler.Enabled
           if systemCfg.Scheduler.force_rank_1_for_approx_wsr_calculations
               num_layers_for_wsr_approx= 1;
           end
           if systemCfg.ChannelSettings.Quadriga_used
               [idx_AS_ue, idx_AS_gnb, num_AS] = cluster_gnb_schedule_UE_qd(...
                   systemCfg.Scheduler,idx_asoc_ue_iter,idx_asoc_gnb,SNR_lin,Nu,Nb,Nr_ob,Nt,...
                    overdim_factor,AS_metric,AS_size_max,O_complexity,...
                    lkg_exponent,max_lkg_allowed, updates_overiters.ue_priority, clustering_type,...
                    coordg_x, coordg_y, systemCfg.gNB.BBU.Cell_Association, ...
                    num_layers_for_wsr_approx, ffinfo_cluster);
           else
               if systemCfg.Traffic_model.Enabled
                   % check each UE's buffer and remove UEs without arrival packets in
                   % current time slot
                   updates_overiters = traffic_check(iter,updates_overiters,systemCfg);
                   idx_asoc_ue = updates_overiters.traffic_model.idx_asoc_ue_bursty;
                   %idx_asoc_ue = idx_asoc_ue(~cellfun('isempty',idx_asoc_ue));
               end
               [idx_AS_ue, idx_AS_gnb, num_AS] = cluster_gnb_schedule_UE(...
                    systemCfg.Scheduler,idx_asoc_ue,idx_asoc_gnb,snrlin_cl,Nu,Nb,Nr_ob,Nt,...
                    overdim_factor,AS_metric,AS_size_max,O_complexity,...
                    lkg_exponent,max_lkg_allowed, updates_overiters.ue_priority, clustering_type,...
                    coordg_x, coordg_y,systemCfg.NetworkLayout, systemCfg.gNB.BBU.Cell_Association, ...
                    num_layers_for_wsr_approx, ffinfo_cluster,Nx,...
                    systemCfg.HetNet.hotspot_numHotspots,systemCfg.HetNet.hotspot_reuse1, iter);
           end
       end
       NpI=getNoiseVector(idx_AS_ue,idx_AS_gnb,SNR_lin).';
       PLs_scale =PLs;
       if systemCfg.ChannelSettings.Quadriga_used
           [SINR_lin,Cap,Rate,BitsPerIter, dn, tm, idx_AS_ue,idx_AS_gnb, ...
                sim_time, sim_time2, Htot_rc,Htot_h_rc, Ht_snr_rc, Ht_snr_h_rc, Indx_ConsRetx_iter, ...
                prevServingAP_iter, reTX_TabuList] = run_one_iter_qd_main_v5(...
                Htot,Htot_h, Ht_snr, Ht_snr_h, iter,...
                Htot_rc,Htot_h_rc, Ht_snr_rc, Ht_snr_h_rc, systemCfg, sim_time, ...
                sim_time2, dn, tm, NpI, PLs_scale, gnbUEBeam_asoc, UEsWithData, ...
                updates_overiters.ue_priority, noise_floor, idx_AS_ue,idx_AS_gnb, ...
                Nr,Nr_ob,L,Nt,Nt2,Nu,Nb,Nfreq,BW,num_layers_for_wsr_approx,...
                gap2cap_lin, isPerantennaSNR, Rx_coloring, chanPerturb, ...
                Indx_ConsRetx, prevServingAP, reTXsuccessiveFailure, reTX_TabuList, ...
                gnbUEBeam_asoc_org_diffQ, p_store, pa_store, Ht_all, Ht_snr_all);
             prevServingAP = [prevServingAP prevServingAP_iter];
             Indx_ConsRetx = [Indx_ConsRetx Indx_ConsRetx_iter];
        else
            % Set some Quadriga specific params to zero
            Nt2= NaN;
            if isUL
                [SINR_lin,Cap, RC_output, updates_overiters]=run_one_iter_UL(systemCfg,...
                    Ht_snr,Ht_snr_h,num_AS,...
                    idx_AS_ue,idx_AS_gnb,updates_overiters);
            else
                [SINR_lin,Cap, ~, inv_cond, ~, inrmacro, inrsc2macro, updates_overiters, RC_output]=run_one_iter(...             %%% Wanlu: SINR_lin = Nu by Nfreq, one value for each (user,tone)
                    Htot,Htot_h,Ht_snr,Ht_snr_h,SNR_cap,...
                    Nr,Nr_ob,L,Nt,Nt2,Nu,Nb,Nfreq,num_AS,...
                    idx_AS_ue,idx_AS_gnb,...
                    chanPerturb,Rx_coloring,systemCfg.gNB.BBU.BeamForming,...
                    gap2cap_lin,isPerantennaSNR, AS_size_max,...
                    updates_overiters.ue_priority, num_layers_for_wsr_approx,NpI, PLs_asoc,systemCfg.UE.Demapper,...
                    systemCfg.gNB.BBU.Cell_Association.coord_silencing_post_clustering, systemCfg.HetNet.hotspot_isTrackMacroInterf,...
                    systemCfg.gNB.BBU.Cell_Association.clustering_init_off_allow_overlap,Nb_gNB, ...
                    systemCfg, updates_overiters);
               %if systemCfg.gNB.Rate_Control.RankSelEnable == 1
               %sinr_offset_ues_iters(:,:,iter) = sinr_offset_ues;
               %else
               %sinr_offset_ues_iters(:,:,iter) = sinr_offset_ues(:,L);
               %end            
               detection_ues_iters(:,iter,d) = detect_outcome;
           end
           if SC_outdoor.isTrackMacroInterf
               INR_macro_iters(:,iter)=inrmacro;
               INR_sc2macro_iters(:,iter)=inrsc2macro;
           end
       end
       %-----------------collect stats per iter----------
       if systemCfg.ChannelSettings.Quadriga_used
           [SINR_tot_lin(:,iter), SINR_tot(:,iter),Cap_tot(:,iter), Tput_tot(:,iter), Bits_tot(:, iter), Rate_tot(:,iter)]=...
               collect_stats_per_iter_qd(SINR_lin,Cap,BitsPerIter,Rate,Nfreq, Sys_Tput_factor,...
               L, systemCfg.UE.rate_per_layer_cap, systemCfg.techStandard.st80211ad);
       else
           [SINR_tot_lin(:,iter), SINR_tot(:,iter),Cap_tot(:,iter), Tput_tot(:,iter), Rate_oh_tot(:,iter), ...
               SINR_lin_eff_tot(:,iter), mcs_ues_tot(:,iter), rank_ues_tot(:,iter), Rate_rec_oh_tot(:,iter), detect_outcome_tot(:,iter)]=...
               collect_stats_per_iter(SINR_lin,Cap,Nfreq, Sys_Tput_factor,...
               L, systemCfg.UE.rate_per_layer_cap, RC_output);
            % [find(SINR_tot_lin(:,iter)>0) Cap_tot(SINR_tot_lin(:,iter)>0,iter)] % for debug
       end
       idx_AS_ue_log{iter} = idx_AS_ue;
       idx_AS_gnb_log{iter} = idx_AS_gnb;
   
   end %  iter loop
    
    %---------------Collect stats per drop-----
    stats=collect_stats_per_drop(stats, SNR_asoc,SNR,PLs_asoc,Prx_asoc,...
        Cap_tot,num_iter,SINR_tot_lin,SINR_tot, Tput_tot, Rate_oh_tot, SINR_lin_eff_tot, Rate_rec_oh_tot, updates_overiters);
    
    stats.mcs{d} = mcs_ues_tot;
    stats.rank{d} = rank_ues_tot;
    sinr_offset_ues_drops{d} = sinr_offset_ues_iters;
    stats.detect{d} = detect_outcome_tot;
    
    stats.LOS=[stats.LOS;LOS];
    if and(Nx~=0,rrue_type=="dim_shp_rept")
        stats.snr_asoc_g2x=[stats.snr_asoc_g2x;SNR_asoc_g2x];
        stats.snr_asoc_x2u=[stats.snr_asoc_x2u;SNR_asoc_x2u];
    end

    if and(Nx~=0,rrue_type=="dim_shp_rept")
        stats.snr_asoc_g2x=[stats.snr_asoc_g2x;SNR_asoc_g2x];
        stats.snr_asoc_x2u=[stats.snr_asoc_x2u;SNR_asoc_x2u];
    end
    
    if SC_outdoor.isTrackMacroInterf
        inter_temp1 = INR_macro_iters(INR_macro_iters~=0);
        stats.macro2sc_inr{d} = 10*log10(10*log10(abs(inter_temp1)));
        inter_temp2 = INR_sc2macro_iters(INR_sc2macro_iters~=0);
        stats.sc2macro_inr{d} = 10*log10(10*log10(abs(inter_temp2)));
    end
    
   if isSave >= 2  % detailed logging
       stats.idx_AS_ue{d}=idx_AS_ue_log;
       stats.idx_AS_gnb{d}=idx_AS_gnb_log;
       stats.cap_iter_drop{d}=Cap_tot;
       if systemCfg.ChannelSettings.Quadriga_used
           stats.BlockedAPs_drop{d} = BlockedAPs;
           stats.Rate_iter_drop{d}=Rate_tot;
           stats.Bits_iter_drop{d}=Bits_tot;
           stats.SINR_iter_drop{d}= SINR_tot;
           stats.traffic{d} = logTrafficModel(tm,d);
           stats.UEassocgNbMTX{d} = idx_asoc_ue;
           stats.ConsectivereTXFailureFlag{d} = Indx_ConsRetx;
           stats.prevServingAP_drop{d} = prevServingAP;
           stats.UEassocTopAPs{d} =gnbUEBeam_asoc_org;
           stats.UEassocTopAPsdiffQ{d} =gnbUEBeam_asoc_org_diffQAll;
       end
   end
    
    % save per drop (useful for long sims)
    saveSimulationData(systemCfg, stats, configureXML, 'temp')
end

% one final save
saveSimulationData(systemCfg, stats, configureXML, 'final')
closeFiles(fids_struct);

