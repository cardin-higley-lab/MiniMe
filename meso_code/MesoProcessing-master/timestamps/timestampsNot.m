function [returnOn, returnOff] = timestampsNot(timestampsOn,timestampsOff,varargin)
%TIMESTAMPNOT returns timetamps corresponding to the inverse of the input
 % Swaps "(true,true,false,true) to (false,false,true,false).
 % Can provide start and end time with optional arguments
 % EXAMPLE
 % [returnOn, returnOff] = ...
 % timestampsNot(timestampsOn,timestampsOff,'StartTime',1,'EndTime',100000);
 returnOn = timestampsOff(1:end-1);
 returnOff = timestampsOn(2:end);
 
 
 StartTimeInd = find(strcmp('StartTime',varargin),1);
 EndTimeInd = find(strcmp('EndTime',varargin),1);
 if ~isempty(StartTimeInd)
    returnOn =[varargin{StartTimeInd+1} returnOn];
    if ~ isempty(timestampsOn)
        returnOff =[timestampsOn(1) returnOff];
    end
 end
 if ~isempty(EndTimeInd)
     if ~ isempty(timestampsOff)
        returnOn =[returnOn timestampsOff(end)];
     end
    returnOff =[returnOff varargin{EndTimeInd+1}];
 end
end