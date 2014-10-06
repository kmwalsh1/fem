function genPointsLoads(InputName, NormName, IsppaNorm, sym, LCID);
% function genPointsLoads(InputName, NormName, IsppaNorm, sym, LCID);
%
% INPUTS:
% InputName (string) - dyna*.mat file to process
% NormName (string) - dyna*.mat file with known Isppa 
% IsppaNorm (float) - Isppa value for the normalization file (W/cm^2)
% PulseDuration (float) - pulse duration (us)
% sym (string) - symmetry boundary conditions
%                'q' -> quarter symmetry
%                'h' -> half symmetry
%                'none' -> no symmetry
% LCID - load curve ID for the point loads
% 
% OUTPUTS:
% PointLoads_a*.dyn is written to the CWD
%
% EXAMPLE:
% genPointsLoads('dyna_ispta_att0.5.mat', 'dyna_ispta_att0.5.mat', ...
%                1000, 'q', 1);
%

%% TODO %%
% [ ] check for uniform element volume (instead of user input): calculate if
%     uniform, o/w assume non-uniform
% [ ] add calls to Python code to generate element volume file

DEBUG = 0;

% node tolerance to search for center line in the lateral
% dimension
LatTol = 1e-3;  % cm

% find the Isppa value that Field II solved for in the
% normalization case (limiting the search to the focal zone,
% +/- 25% the focal depth)
AxialSearch = 0.25; % percentage of the focal depth to search
										%for the Isppa value
load(NormName);
NormIntensity = intensity;
mpn = FIELD_PARAMS.mpn;
NormFocalDepth = FIELD_PARAMS.focus(3)*100;  % convert m -> cm

% check to make sure nodes exist at lat == 0 for the push
if(isempty(find(abs(mpn(:,3)) < LatTol))),
    keyboard
    error('lat = 0 nodes missing for rad force excitation');
end;

[NormAx, NormFZ] = extractAxisIntensity(NormName, 3, [0 0]);

% what is the Isppa value that field has solved
NormFieldIsppa = max(NormFZ)

% make plot of the intensity profile to make sure that
% everything makes sense
if DEBUG,
    figure;
    plot(abs(mpn(NormFZ,4)),NormIntensity(NormFZ),'-kx');
    hold on;
end;

% find normalization max in desired alpha
load(InputName);
InputIntensity = intensity;
mpn = FIELD_PARAMS.measurementPointsandNodes;
FocalDepth = FIELD_PARAMS.focus(3)*100;  % convert m -> cm

[Ax, FZ] = extractAxisIntensity(InputName, 3, [0 0]);

% what is the Isppa value that field has solved
FieldIsppa = max(FZ)

if DEBUG,
    % add this one to the plot
    plot(abs(mpn(FZ,4)),InputIntensity(FZ),'-rx');
    xlabel('Depth (cm)');
    ylabel('Field Intensity');
    title('Comparison of Field Intensity Profiles');
    legend('Normalization','Input','Location','Best');
    legend boxoff;
end;

% normalize InputIntensity
InputIntensity = InputIntensity./FieldIsppa;

% toss intensities below 5% of Isppa
InputIntensity=InputIntensity.*(InputIntensity>=0.05);

% now zero out values near the transducer face b/c they
% violated the farfield assumption in field
z0=find(abs(mpn(:,4)) < 0.001);
InputIntensity(z0)=0;

% scale the Field intensities relative to the known intensity value for the
% normalization data
Field_I_Ratio = FieldIsppa/NormFieldIsppa
ScaledIntensity = InputIntensity*IsppaNorm*Field_I_Ratio;

SoundSpeed = FIELD_PARAMS.soundSpeed*100;  % convert m/s -> cm/s

% open up an ASCII file for writing point loads and initial temperatures and
% add comments telling the user what runtime parameters were used
RunTimeDate = date;

fout_filename = sprintf('PointLoads-f%.2f-F%.1f-FD%.3f-a%.2f.dyn', FIELD_PARAMS.Frequency, FIELD_PARAMS.Fnum, ...
                                                                   FIELD_PARAMS.focus(3), FIELD_PARAMS.alpha);
foutload = fopen(fout_filename,'w');
fprintf(foutload,'*LOAD_NODE_POINT\n');
fprintf(foutload,'$ Generated by genPointLoads on %s\n',RunTimeDate);
fprintf(foutload,'$ Normalization File:\n');
fprintf(foutload,'$ %s\n',NormName);
fprintf(foutload,'$ Normalization Isppa = %.1f W/cm^2\n',IsppaNorm);

% convert alpha -> Np/cm
AlphaNp = FIELD_PARAMS.alpha*FIELD_PARAMS.Frequency/8.616;

% solve for point loads
MaxLoad = 0;
MaxBodyForce = 0;
for x=1:length(mpn),
    xcoord = mpn(x,2);
    ycoord = mpn(x,3);
    zcoord = mpn(x,4);
    if (ScaledIntensity(x)~=0 & zcoord~=0) 
        if(isnan(ScaledIntensity(x)) || (ScaledIntensity(x) > (10*IsppaNorm))),
            warning('Excessive intensities are being computed; not writing to output file');
        else,
            NodeID=mpn(x,1);
            % 1 W = 10,000,000 g cm^2/s^2
            ScaledIntensityLoad = ScaledIntensity(x) * 10000000;
            BodyForce = (2*AlphaNp*ScaledIntensityLoad)/SoundSpeed;
            if BodyForce > MaxBodyForce,
                MaxBodyForce = BodyForce;
                ScaledIsppa = ScaledIntensity(x);
                %%%%%%%%%%%%%%%%%% CHANGE BASED ON MESH UNIFORMITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isa(ElementVolume,'float'),
                PointLoad = BodyForce * ElementVolume;
            elseif isa(ElementVolume,'char'),
                PointLoad = BodyForce * NodeVolumes(find(NodeVolumes(:,1) == NodeID),2);
            end;
                %%%%%%%%%%%%%%%%%% CHANGE BASED ON MESH UNIFORMITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch sym,
                case 'q'
                    % if the load is on the symmetry axis (x = y = 0), then divide by 4; if
                    % not, check if it is on a symmetry face (x = 0 || y = 0), then divide
                    % by 2
                    if(abs(xcoord) < 1e-4 && abs(ycoord) < 1e-4),
                        PointLoad = PointLoad / 4;
                    elseif(abs(xcoord) < 1e-4 || abs(ycoord) < 1e-4),
                        PointLoad = PointLoad / 2;
                    end;
                case 'h'
                    % if the load is on the symmetry face (x=0), then divide by 2
                    if(abs(xcoord) < 1e-4),
                        PointLoad = PointLoad / 2;
                    end;
                otherwise
                    disp('No symmetry load scaling performed.')
            end;
            if(abs(PointLoad) > MaxLoad),
                MaxLoad = abs(PointLoad);
            end;
            % write point load data (negative to point in -z direction in the
            % dyna model)
            fprintf(foutload,'%d,3,%d,%0.2e,0 \n', NodeID, LCID, -PointLoad);
        end;
    end;
end;

disp(sprintf('Point loads are pointed in the -z direction.\n'));

disp(sprintf('Isppa = %.1f W/cm^2\n',ScaledIsppa));
disp(sprintf('MaxLoad = %.4f g cm/s^2\n',MaxLoad));

fprintf(foutload,'*END');
fclose(foutload);

disp('The following output file has been written:');
disp(fout_filename);

disp(sprintf('This can be included in the LS-DYNA input decks.\n'));
