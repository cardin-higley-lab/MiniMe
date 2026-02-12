function spike2_data = just_timestamps(dataSmrxFile, cedpath, vis_isi)
%% extract timestamps from smrx file and align imaging data to spike2 data
    disp('Extracting smrx timestamps');
    spike2_clayton = loadSpike2Data(dataSmrxFile, cedpath);
    disp('processing smrx timestamps');
    spike2_data = processSpike2Data(spike2_clayton, vis_isi);
end