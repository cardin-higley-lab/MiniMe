function fig = stateplots(spike2_data, pc, stateOn, stateOff, outputPath, figName, fig_title)

fig = figure();
ax(1) = subplot(211);
plot(spike2_data.analog_signal_time_vect, spike2_data.wheelSpeed);
ylabel('wheel speed (cm/s)');
xline(stateOn, 'g');
xline(stateOff, 'r');
title(sprintf('%s', fig_title));

ax(2) = subplot(212);
maxlength = min(length(spike2_data.pupilFrameOnTimestamps), length(pc));
plot(spike2_data.pupilFrameOnTimestamps(1:maxlength), pc(1:maxlength));
ylabel('facemap pc1');
xline(stateOn, 'g');
xline(stateOff, 'r');
linkaxes(ax,'x');


saveas(fig, fullfile(outputPath, 'stateplots', figName));

end