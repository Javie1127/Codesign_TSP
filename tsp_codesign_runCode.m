function tsp_codesign_runCode(configureXML)

%read XML config file
if nargin==0
    clear;
    configureXML = 'Example_1/Case_1/system_cfg.xml';
end
tmpStruct = xml2struct(configureXML);
systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);
rng(systemCfg.System_Parameters.RAND_SEED);
num_iter = systemCfg.System_Parameters.NUM_DROPS;
% Initialization before iterations
[fdcomm, radar, radar_comm,cov] = tsp_initializations(systemCfg);
%stats = init_stats();s
%%% 
for d=1:num_iter   
    %% BCD-AP algorithm 
    [fdcomm_op, radar_op] = tsp_BCD_AP(fdcomm, radar, radar_comm, cov);
    
end
save('case_1.mat','fdcomm_op','radar_op');
% 
% % one final save
% saveSimulationData(systemCfg, stats, configureXML, 'final')
% closeFiles(fids_struct);
end

