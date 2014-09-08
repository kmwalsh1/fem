function [impulseResponse]=formatExpImpResp(FIELD_PARAMS)
% function [impulseResponse]=formatExpImpResp(FIELD_PARAMS)
%
% Format the raw oscilloscope data to create an impulse response vector that
% can be used by Field II
%
% INPUT:
%   FIELD_PARAMS (struct) - structure of run time parameters
%
% OUTPUT:
%   impulseResponse (float vector) - experimentally-measured impulse response
%   that has been centered and sampled correctly
%
% Mark 06/20/07

SampFreq = FIELD_PARAMS.samplingFrequency;

% read in the raw data from the oscilloscope
if(~isfield(FIELD_PARAMS.probeStruct.impulse_response, 'time') ||...
   ~isfield(FIELD_PARAMS.probeStruct.impulse_response, 'voltage'))
    error(sprintf('The experimentally measured impulse response doesn''t exist in %s.json.\n',FIELD_PARAMS.Transducer));
end;

TimePulse = FIELD_PARAMS.probeStruct.impulse_response.time;
% normalize the pulse data
Pulse = FIELD_PARAMS.probeStruct.impulse_response.voltage...
        ./max(FIELD_PARAMS.probeStruct.impulse_response.voltage);
% center the time axis around the max intensity
[MaxPulse,MaxPulseIndex]=max(Pulse);
TimePulse = TimePulse - TimePulse(MaxPulseIndex);

% re-sample the data to match the Field II sampling frequency
NewTime = min(TimePulse):1/SampFreq:max(TimePulse);
NewPulse = interp1(TimePulse,Pulse,NewTime);

% find the indices of NewTime to use for the Field II impulse response
StartIndex = min(find(NewTime > -4e-7));
StopIndex = min(find(NewTime > 4e-7));

impulseResponse = NewPulse(StartIndex:StopIndex);
