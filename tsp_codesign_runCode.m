function tsp_codesign_runCode(configureXML)

%read XML config file
if nargin==0
    clear;
    configureXML = 'Example_1/system_cfg.xml';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
rng(systemCfg.System_Parameters.RAND_SEED);
num_iter = systemCfg.System_Parameters.NUM_DROPS;
% Initialization before iterations
[fdcomm, radar, radar_comm,cov] = tsp_initializations(systemCfg);

%%% 
for d=1:num_iter
    % iniializing channel models
    [fdcomm, radar, radar_comm] = tsp_channel_inis(fdcomm,radar,radar_comm);
    [fdcomm, radar, cov] = tsp_precoder_ini(radar,fdcomm,radar_comm,cov);
    % save per drop (useful for long sims)
    saveSimulationData(systemCfg, stats, configureXML, 'temp')
end

% one final save
saveSimulationData(systemCfg, stats, configureXML, 'final')
closeFiles(fids_struct);
end

