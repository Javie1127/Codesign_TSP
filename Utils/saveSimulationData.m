function saveSimulationData(systemCfg, stats, fdcomm,radar, configureXML, saveType)
isSave = systemCfg.System_Parameters.isSave;
if isSave
    [filepath,name,~] = fileparts(configureXML);
    precoder_type = systemCfg.Cases.CodingSchemes.Comm_Precoding_Scheme;
    coding_type = systemCfg.Cases.CodingSchemes.Radar_Coding_Scheme;
    Initializations = systemCfg.Cases.Initializations;
    if strcmpi(saveType,'final')
        mat_filename = sprintf('sim_output_%s_%s_%s_%s.mat', precoder_type,coding_type,...
    Initializations,systemCfg.iteration_case);
    else
        mat_filename = sprintf('temp_sim_output_%s_%s_%s_%s_%d.mat', precoder_type,coding_type,...
    Initializations,systemCfg.iteration_case,systemCfg.iter_variable);
    end
    savefile = fullfile(filepath,mat_filename);
    save(savefile, 'systemCfg', 'stats','fdcomm','radar');
    fprintf('Simulation output saved to %s\n', savefile);
else
    fprint('No simulation data files were saved: isSave=%d', isSave);
end
end