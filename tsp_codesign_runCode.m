function tsp_codesign_runCode(configureXML)

if isunix
    addpath(genpath([pwd,'/Utils']));
    updateMatlabPath()
else
    addpath(genpath([pwd, '\Utils']));
    updateMatlabPath()
end

%read XML config file
if nargin==0
    clear;
    configureXML = 'Example_1/system_cfg.xml';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
rng(systemCfg.SLS_Top_Cfg.RAND_SEED);
% Initialization before iterations
[fdcomm, radar, radar_comm] = tsp_initializations(systemCfg);

%%% 
for d=1:num_iter
    % iniializing channel models
    [fdcomm, radar, radar_comm] = tsp_channel_inis(fdcomm,radar,radar_comm);
    % save per drop (useful for long sims)
    saveSimulationData(systemCfg, stats, configureXML, 'temp')
end

% one final save
saveSimulationData(systemCfg, stats, configureXML, 'final')
closeFiles(fids_struct);
end

