function tsp_codesign_runCode(configureXML)

addpath(genpath([pwd, '\Utils']));
updateMatlabPath()


%read XML config file
if nargin==0
    clear;
    configureXML = 'Example_1\EW_EI_SNR_50_dB\system_cfg.xml';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
rng(systemCfg.System_Parameters.RAND_SEED);
num_iter = systemCfg.System_Parameters.NUM_Channel_REALIZATIONS;
num_rept = systemCfg.System_Parameters.NUM_DROPS;
% Initialization before iterations
[fdcomm_ini, radar_ini, radar_comm_ini,cov_ini] = tsp_initializations(systemCfg);
fprintf('current case: %s_%s_%s_%s_%d_%d ',fdcomm_ini.precoder_type,radar_ini.coding_type,...
    fdcomm_ini.Initializations,fdcomm_ini.Weights,fdcomm_ini.DL_num,fdcomm_ini.UL_num);


I_total = zeros(radar_ini.ell_max,num_iter);
I_total_UL = zeros(radar_ini.ell_max,num_iter);
I_total_DL = zeros(radar_ini.ell_max,num_iter);
I_total_radar = zeros(radar_ini.ell_max,num_iter);
precoder_type = fdcomm_ini.precoder_type;
radar_coding_type = radar_ini.coding_type;
parfor d=1:num_iter % num_iter channel realizations
    %% initializing channel models
    [fdcomm_d, radar_d, radar_comm_d] = tsp_channel_inis(fdcomm_ini,radar_ini,radar_comm_ini);
    %stats_ini = init_stats();
    I_tmp_track = 0;
    fdcomm_op_d = fdcomm_d;
    radar_op_d = radar_d;
    for q = 1:num_rept % num_rep > 1 for the random initialization only otherwise num_rep = 1
        disp(d);
        switch precoder_type
            case 'Proposed'
                fdcomm = tsp_precoder_ini(fdcomm_d);
            case 'BD'% DL only
                fdcomm = tsp_codesign_BD_DL(fdcomm_d);
            case 'NSP'% DL only UL uniform precoding
                fdcomm = tsp_codesign_NSP(fdcomm_d,radar_comm_d);
            case 'Uniform' % UL Only
                fdcomm = tsp_codesign_Uniform_Precoding_UL(fdcomm_d);
        end
        switch radar_coding_type
            case 'Proposed'
                radar=tsp_radar_code_ini(radar_d);
            case 'Random'
                radar = randomcode(radar_d);
            case 'All-one'
                radar = radar_all_one(radar_d);
        end
        %% Initializing the covariance matrices
        cov = tsp_covmat_rev(fdcomm,radar,radar_comm_d,cov_ini);
        %% Initializing the MMSE matrices
        fdcomm = tsp_Comm_MMSE_rev(fdcomm,radar,cov);
        radar = tsp_radar_MMSE_rev(radar,cov);
        %% Initializing the performance measures
        radar = Xi_radar(radar);
        K = radar.codelength;
        for k = 1:K
            fdcomm = Xi_comm_k(fdcomm,k);
        end
        %% BCD-AP algorithm 
        [fdcomm_tmp, radar_tmp] = tsp_BCD_AP(fdcomm, radar, radar_comm_d, cov);
        %% Update Stats
        %stats_temp = collect_stats(fdcomm_op_d,radar_op_d,stats_ini);
        if max(fdcomm_tmp.I_total_op)>I_tmp_track
            fdcomm_op_d = fdcomm_tmp;
            radar_op_d = radar_tmp;
            I_tmp_track = max(fdcomm_tmp.I_total_op);
        end
    end
    saveSimulationData_parfor(systemCfg, fdcomm_op_d,radar_op_d, configureXML,'tmp');
    I_total(:,d) = fdcomm_op_d.I_total_op;
    I_total_DL(:,d) = fdcomm_op_d.I_DL_op;
    I_total_UL(:,d) = fdcomm_op_d.I_UL_op;
    I_total_radar(:,d) = radar_op_d.I_radar_op;
end
%stats_final = stats_temp;
stats = collect_stats_final_parfor(I_total,I_total_DL,I_total_UL,I_total_radar);
saveSimulationData(systemCfg,stats, fdcomm_ini,radar_ini, configureXML,'final');
% 
% % one final save
% 
% closeFiles(fids_struct);
end

