function tsp_codesign_runCode_multi_case_v3(configureXML)

addpath(genpath([pwd, '\Utils']));
updateMatlabPath()


%read XML config file
if nargin==0
    clear;
    configureXML = 'Example_4\Proposed_Proposed_Co\system_cfg.xml';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
rng(systemCfg.System_Parameters.RAND_SEED);
num_iter = systemCfg.System_Parameters.NUM_Channel_REALIZATIONS;


switch systemCfg.iteration_case
    case 'SNR_radar'
        iter_variables = systemCfg.Powers.SNR_radar;
    case 'CNR'
        iter_variables = systemCfg.Powers.CNR;
    case 'UL_SNR'
        iter_variables = systemCfg.Powers.SNR_UL;
    case 'DL_SNR'
        iter_variables = systemCfg.Powers.SNR_DL;
    case 'DL_UE'
        iter_variables = systemCfg.Cases.FD.Num_DL_UE;
    case 'UL_UE'
        iter_variables = systemCfg.Cases.FD.Num_UL_UE;
    case 'PAR'
        iter_variables = systemCfg.Powers.radar_PAR;
end
stats_total = cell(length(iter_variables),1);
for ii = 1:length(iter_variables)
    systemCfg_ii = systemCfg;
    if strcmp(systemCfg_ii.iteration_case, 'SNR_radar')
        systemCfg_ii.Powers.SNR_radar = iter_variables(ii) ;
    elseif strcmp(systemCfg_ii.iteration_case, 'CNR')
        systemCfg_ii.Powers.CNR = iter_variables(ii) ;
    elseif strcmp(systemCfg_ii.iteration_case, 'UL_SNR')
        systemCfg_ii.Powers.SNR_UL = iter_variables(ii) ;
    elseif strcmp(systemCfg_ii.iteration_case, 'DL_SNR')
        systemCfg_ii.Powers.SNR_DL = iter_variables(ii) ;
    elseif strcmp(systemCfg_ii.iteration_case, 'DL_UE')
        systemCfg_ii.Cases.FD.Num_DL_UE=iter_variables(ii);
    elseif strcmp(systemCfg_ii.iteration_case, 'UL_UE')
        systemCfg_ii.Cases.FD.Num_UL_UE=iter_variables(ii);
    elseif strcmp(systemCfg_ii.iteration_case, 'PAR')
        systemCfg_ii.Powers.radar_PAR = iter_variables(ii) ;
    end
    systemCfg_ii.iter_variable = iter_variables(ii);
    % Initialization before iterationssystemCfg_ii.Powers.SNR_radar = iter_variable(ii) ;
    [fdcomm_ini, radar_ini, radar_comm_ini,cov_ini] = tsp_initializations(systemCfg_ii);
    I_total = zeros(radar_ini.ell_max,num_iter);
    I_total_UL = zeros(radar_ini.ell_max,num_iter);
    I_total_DL = zeros(radar_ini.ell_max,num_iter);
    I_total_radar = zeros(radar_ini.ell_max,num_iter);
    precoder_type = fdcomm_ini.precoder_type;
    radar_ii = cell(num_iter,1);
    fdcomm_ii = cell(num_iter,1);
    parfor d=1:num_iter % num_iter channel realizations
        %% initializing channel models
        [fdcomm_d, radar_d, radar_comm_d] = tsp_channel_inis(fdcomm_ini,radar_ini,radar_comm_ini);
        %stats_ini = init_stats();
        disp(d);
        radar_coding_type = radar_ini.coding_type;
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
                %radar = randomcode(radar_d);
            case 'Random'
                radar = randomcode(radar_d);
                radar = tsp_Nearest_PAR(radar);
            case 'All-one'
                radar = radar_all_one(radar_d);
        end
        %% Initializing the covariance matrices
        cov = tsp_covmat_rev(fdcomm,radar,radar_comm_d,cov_ini);
        %% Initializing the WMMSE Rxs
        fdcomm = tsp_Comm_MMSE_rev(fdcomm,radar,cov);
        radar = tsp_radar_MMSE_rev(radar,cov);
%             %% Initializing the performance measures
%             radar = Xi_radar(radar);
%             K = radar.codelength;
%             for k = 1:K
%                 fdcomm = Xi_comm_k(fdcomm,k);
%             end
        %% BCD-AP algorithm 
        [fdcomm_op_d, radar_op_d] = tsp_BCD_AP(fdcomm, radar, radar_comm_d, cov);
        I_total(:,d) = fdcomm_op_d.I_total_op;
        I_total_DL(:,d) = fdcomm_op_d.I_DL_op;
        I_total_UL(:,d) = fdcomm_op_d.I_UL_op;
        I_total_radar(:,d) = radar_op_d.I_radar_op;
        fdcomm_ii{d,1} = fdcomm_op_d;
        radar_ii{d,1} = radar_op_d;
    end
    stats = collect_stats_final_parfor(I_total,I_total_DL,I_total_UL,I_total_radar);
    saveSimulationData(systemCfg_ii,stats,fdcomm_ii,radar_ii, configureXML,'tmp');
    stats_total{ii} =  stats;
end
saveSimulationData(systemCfg,stats_total, fdcomm_ini,radar_ini, configureXML,'final');

% 
% % one final save
% 
% closeFiles(fids_struct);
end

