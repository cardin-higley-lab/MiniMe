function [outputArg1,outputArg2] = locoFaceConjunctiveStates(timestamp1On,timestamp1Off,timestamp2On,timestamp2Off)
%locoFaceConjunctiveStates returns the timestamps corresponding to the
%of (high face and low locomotion), (high or low face and locomotion), and
%(low face and low locomotion)
%   Detailed explanation goes here



for i = 1:length(timestamp1Off)
    for ii = 1:length(timestamp2Off)
        [andOn,andOff] = timestampAnd(timestamp1On(i),timestamp1Off(i),timestamp2On(ii),timestamp2Off(ii));
        if ~isempty(andOn)
            conjunctiveTimestamps.quiescenceOn
    end
end

