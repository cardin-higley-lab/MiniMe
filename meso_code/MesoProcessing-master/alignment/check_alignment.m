figure
ax(1) = subplot(311);
plot(spike2_data.analog_signal_time_vect,spike2_data.wheelSpeed); 
xlabel('wheel speed');


ax(2) = subplot(312);
l = min(length(mesoTimestamps)/2, length(dFoF_parcels));
plot(mesoTimestamps(1:2:l), nanmean(dFoF_parcels(:,1:l)',2));
xlabel('mean parcel trace - mesoTimestamps');


ax(3) = subplot(313);
l = min(length(spike2_data.blueOnTimestamps), length(dFoF_parcels));
plot(spike2_data.blueOnTimestamps(1:l), nanmean(dFoF_parcels(:,1:l)',2));
xlabel('mean parcel trace - blueOnTimestamps');

linkaxes(ax,'x')