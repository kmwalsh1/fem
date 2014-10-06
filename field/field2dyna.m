function []=field2dyna(NodeName,alpha,Fnum,focus,Frequency,Transducer,Impulse)
% function []=field2dyna(NodeName,alpha,Fnum,focus,Frequency,Transducer,Impulse)
%
% INPUT:
% NodeName (string) - file name to read nodes from (e.g., nodes.dyn); needs to
% be a comma-delimited file with header/footer lines that
% start with *
% alpha - 0.5, 1.0, etc. (dB/cm/MHz)
% Fnum - F/# (e.g. 1.3)
% focus - [x y z] (m) "Field" coordinates
% Frequency - excitation frequency (MHz)
% Transducer (string) - 'vf105','vf73'
% Impulse (string) - 'gaussian','exp'
%
% OUTPUT: dyna_ispta*.mat file is saved to CWD
%
% EXAMPLES:
% field2dyna('nodes.dyn',0.5,1.3,[0 0 0.02],7.2,'vf105','gaussian');
%

mpn = read_mpn(NodeName);

% check to see if there are nodes in the x = y = 0 plane to make sure that
% "on axis" intensities are properly captured
check_on_axis(mpn(:,2:4));

% invert the z axis
mpn(:,4)=-mpn(:,4);

% switch x and y so plane of symmetry is elevation, not lateral
mpn(:,2:3)=[mpn(:,3) mpn(:,2)];
  
% convert from cm -> m
mpn(:,2:4)=mpn(:,2:4)/100;

% create a variable structure to pass to dynaField
FIELD_PARAMS.mpn = mpn;
FIELD_PARAMS.alpha = alpha;
FIELD_PARAMS.Fnum = Fnum;
FIELD_PARAMS.focus = focus;
FIELD_PARAMS.Frequency = Frequency;
FIELD_PARAMS.Transducer = Transducer;
FIELD_PARAMS.Impulse = Impulse;

% below are hard-coded constants (transducer independent)
FIELD_PARAMS.soundSpeed=1540;
FIELD_PARAMS.samplingFrequency = 200e6;

% figure out where this function exists to link probes submod
functionDir = fileparts(which(mfilename));
probesPath = fullfile(functionDir, '../probes/fem');
check_add_probes(probesPath);

% check that Field II is in the Matlab search path, and initialize
check_start_Field_II;

% define transducer-independent parameters
set_field('c', FIELD_PARAMS.soundSpeed);
set_field('fs', FIELD_PARAMS.samplingFrequency);

% define transducer-dependent parameters
eval(sprintf('[Th,impulseResponse] = %s(FIELD_PARAMS);', ...
     FIELD_PARAMS.Transducer));

% check specs of the defined transducer
FIELD_PARAMS.Th_data = xdc_get(Th, 'rect');

% figure out the axial shift (m) that will need to be applied to the
% scatterers to accomodate the mathematical element shift due to the lens
FIELD_PARAMS.lens_correction_m = correct_axial_lens(FIELD_PARAMS.Th_data);

% define the impulse response
xdc_impulse(Th, impulseResponse);

% define the excitation pulse
exciteFreq=FIELD_PARAMS.Frequency*1e6; % Hz
ncyc=50;
texcite=0:1/FIELD_PARAMS.samplingFrequency:ncyc/exciteFreq;
excitationPulse=sin(2*pi*exciteFreq*texcite);
xdc_excitation(Th, excitationPulse);

% set attenuation
Freq_att=FIELD_PARAMS.alpha*100/1e6; % FIELD_PARAMS in dB/cm/MHz
att_f0=exciteFreq;
att=Freq_att*att_f0; % set non freq. dep. to be centered here
set_field('att', att);
set_field('Freq_att', Freq_att);
set_field('att_f0', att_f0);
set_field('use_att', 1);
 
% compute Ispta at each location for a single tx pulse
% optimizing by computing only relevant nodes... will assume others are zero
StartTime = fix(clock);
% disp(sprintf('Start Time: %i:%2.0i',StartTime(4),StartTime(5)));
Time = datestr(StartTime, 'HH:MM:SS PM');
disp(sprintf('Start Time: %s', Time))
tic;

numNodes = size(FIELD_PARAMS.mpn, 1);
progressPoints = 0:10000:numNodes;
for i=1:numNodes,
    if ~isempty(intersect(i, progressPoints)),
        disp(sprintf('Processed %.1f%%', i * 100 / numNodes));
    end
    % include the lens correction (axial shift)
    [pressure, startTime] = calc_hp(Th, FIELD_PARAMS.mpn(i,2:4)+FIELD_PARAMS.lens_correction_m);
    intensity(i) = sum(pressure.*pressure);
end

CalcTime = toc; % s
ActualRunTime = CalcTime/60; % min
disp(sprintf('Run Time = %.3f m\n', ActualRunTime));

field_end;

% save intensity file
save(sprintf('dyna-I-f%.2f-F%.1f-FD%.3f-a%.2f.mat', Frequency, Fnum, ...
             focus(3), alpha), 'intensity', 'FIELD_PARAMS');

disp('The next step is to run genPointLoads.');
