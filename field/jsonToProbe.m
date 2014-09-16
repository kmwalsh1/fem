function varargout = jsonToProbe(FIELD_PARAMS);
%function [Th, impulseResponse] = jsonToProbe(FIELD_PARAMS);
%
% Create 'Th' transducer handle and define the impulse response and 
% fractional bandwidth for the given probe for use by dynaField.m;
% impulseResponse will be a structure containing the fractional bandwidth
% and center frequency for acunav probes.

% should change to local fem path, but would need to compile mex file
% before tojson and fromjson functions can be used...
addpath /luscinia/nl91/matlab/matlab-json/

eval(sprintf('probe = jsonToString(''%s.json'');', FIELD_PARAMS.Transducer))

% convert json to workspace variables
for var=fieldnames(probe)'
    eval(sprintf('%s = probe.%s;', char(var), char(var)))
end

% define number of elements
no_elements = (FIELD_PARAMS.focus(3)/FIELD_PARAMS.Fnum)/pitch;
no_elements = floor(no_elements);

if (exist('probe.no_elements', 'var') && no_elements > probe.no_elements)
    no_elements = probe.no_elements;
end

% the function fromjson converts json arrays into matlab cell arrays, so
% here we convert them into a numeric array.
if (isfield(probe.impulse_response, 'time') && isfield(probe.impulse_response, 'voltage'))
    probe.impulse_response.time = cell2mat(probe.impulse_response.time);
    probe.impulse_response.voltage = cell2mat(probe.impulse_response.voltage);
end

%% run extra MATLAB commands
commands = fieldnames(FIELD_PARAMS.commands);
for command=commands'
    if ~strcmp(command, 'impulseResponse') && ~strcmp(command, 'Th')
        eval(fprintf('probe.commands.%s', command))
    end
end
FIELD_PARAMS.probeStruct = probe;
% transducer handle and impulse response
fprintf('%s\n', probe.commands.Th)
eval(probe.commands.Th) % define Th

%% setting return arguments
if (isfield(probe.commands, 'impulseResponse'))
    % return Th, impulseResponse
    fprintf('%s\n', probe.commands.impulseResponse)
    eval(probe.commands.impulseResponse) % define impulseResponse
    varargout = cell(1, 2);
    varargout{1} = Th;
    varargout{2} = impulseResponse;
else
    % Th, fractionalBandwidth, centerFrequency
    varargout = cell(1, 3);
    varargout{1} = Th;
    varargout{2} = fractionalBandwidth;
    varargout{3} = centerFrequency;
end

end


function [probe] = jsonToString(probe_file)
%% jsonToString
%  returns the contents of the json file probe_file as a string. fromjson
%  can then be called on the returned string to convert it into a MATLAB
%  structure
eval(sprintf('fid = fopen(''/luscinia/nl91/matlab/probes/json/%s'');',...
     probe_file));
    
probejson = '';
line = fgetl(fid);

while ischar(line)
    probejson = strcat(probejson, line);
    line = fgetl(fid);
end
    
probe = fromjson(probejson);

end
