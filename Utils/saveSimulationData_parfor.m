function saveSimulationData_parfor(systemCfg, fdcomm,radar, configureXML, saveType)
isSave = systemCfg.System_Parameters.isSave;
if isSave
    [filepath,name,~] = fileparts(configureXML);
    iteration_variab
    if strcmpi(saveType,'final')
        mat_filename = sprintf('sim_output_%s_%s_%s_%s_%d_%d_%s_%d.mat', fdcomm.precoder_type,radar.coding_type,...
    fdcomm.Initializations,fdcomm.Weights,fdcomm.DL_num,fdcomm.UL_num,systemCfg.iteration_case);
    else
        mat_filename = sprintf('temp_sim_output_%s_%s_%s_%s_%d_%s_%d.mat', fdcomm.precoder_type,radar.coding_type,...
    fdcomm.Initializations,fdcomm.Weights,fdcomm.DL_num,fdcomm.UL_num,systemCfg.iteration_case,systemCfg.iter_variable);
    end
    savefile = fullfile(filepath,mat_filename);
    save(savefile, 'systemCfg', 'fdcomm','radar');
    fprintf('Simulation output saved to %s\n', savefile);
else
    fprint('No simulation data files were saved: isSave=%d', isSave);
end
end