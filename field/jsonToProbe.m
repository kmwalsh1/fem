function [Th, impulseResponse] = jsonToProbe(FIELD_PARAMS);
%function [Th, impulseResponse] = jsonToProbe(FIELD_PARAMS);
%
% Create 'Th' transducer handle and define the impulse response and 
% fractional bandwidth for the given probe for use by dynaField.m;

% should change to local fem path, but would need to compile mex file
% before tojson and fromjson functions can be used...
addpath /luscinia/nl91/matlab/matlab-json/

eval(sprintf('probe = jsonToString(''%s.json'');', FIELD_PARAMS.Transducer))
Th = probe;
% need to change formatExpImpResp to take in json formatted time and
% voltage data for the impulse response, rather than the .pul file
end


function [probe] = jsonToString(probe_file)
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