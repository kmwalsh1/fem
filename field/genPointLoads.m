function genPointsLoads(InputName, NormName, IsppaNorm, sym, LCID);
% function genPointsLoads(InputName, NormName, IsppaNorm, sym, LCID);
%
% INPUTS:
% InputName (string) - dyna*.mat file to process
% NormName (string) - dyna*.mat file with known Isppa 
% IsppaNorm (float) - Isppa value for the normalization file (W/cm^2)
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

% find the Isppa value that Field II solved for in the normalization case 
[NormIntensity, NormFocalDepth, mpn] = loadDynaMat(NormName);

checkLat0(mpn);

[NormAx, NormFZ] = extractAxisIntensity(NormName, 3, [0 0]);

% what is the Isppa value that field has solved
NormFieldIsppa = max(NormFZ)

% find normalization max in desired alpha
[InputIntensity, FocalDepth, mpn] = loadDynaMat(InputName);

[Ax, FZ] = extractAxisIntensity(InputName, 3, [0 0]);

% what is the Isppa value that field has solved
FieldIsppa = max(FZ)

% normalize InputIntensity
InputIntensity = InputIntensity./FieldIsppa;

% toss intensities below 5% of Isppa
InputIntensity=InputIntensity.*(InputIntensity>=0.05);

% exclude nodes near xdcr face that could violate far-field assumption
[InputIntensity] = excludeFaceNodes(InputIntensity, mpn);

% scale the Field intensities relative to the known intensity value for the
% normalization data
Field_I_Ratio = FieldIsppa/NormFieldIsppa
ScaledIntensity = InputIntensity*IsppaNorm*Field_I_Ratio;

SoundSpeed = FIELD_PARAMS.soundSpeed*100;  % convert m/s -> cm/s
AlphaNp = FIELD_PARAMS.alpha*FIELD_PARAMS.Frequency/8.616; % convert -> Np/cm

% open up an ASCII file for writing point loads and initial temperatures and
% add comments telling the user what runtime parameters were used

fout_filename = sprintf('PointLoads-f%.2f-F%.1f-FD%.3f-a%.2f.dyn', FIELD_PARAMS.Frequency, FIELD_PARAMS.Fnum, ...
                                                                   FIELD_PARAMS.focus(3), FIELD_PARAMS.alpha);
foutload = fopen(fout_filename,'w');
writeHeader(foutload, NormName, IsppaNorm);


% deterine if the mesh node spacing is uniform
isUniform = checkUniform(mpn(:,2:4));

if isUniform,
    % CALCULATE ELEMENT VOLUME HERE!
else,
    % GENERATE ELEMENT VOLUME FILE
end

% solve for point loads
MaxLoad = 0;
MaxBodyForce = 0;
for x=1:length(mpn),
    xcoord = mpn(x,2);
    ycoord = mpn(x,3);
    zcoord = mpn(x,4);
    if (ScaledIntensity(x)~=0 & zcoord~=0) 
        if(isnan(ScaledIntensity(x)) || (ScaledIntensity(x) > (10*IsppaNorm))),
            warning('Excessive intensities being computed; not writing to output file');
        else,
            NodeID = mpn(x,1);
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
                    % if the load is on the symmetry axis (x = y = 0), then
                    % divide by 4; if not, check if it is on a symmetry face 
                    % (x = 0 || y = 0), then divide by 2
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

if DEBUG,
    % make plot of the intensity profile to make sure that everything makes
    % sense
    figure;
    plot(abs(mpn(NormFZ,4)), NormIntensity(NormFZ),'-kx');
    hold on;
    plot(abs(mpn(FZ,4)), InputIntensity(FZ),'-rx');
    xlabel('Depth (cm)');
    ylabel('Field Intensity');
    title('Comparison of Field Intensity Profiles');
    legend('Normalization','Input','Location','Best');
    legend boxoff;
end;


end

function writeHeader(foutload, NormName, IsppaNorm)
    fprintf(foutload,'*LOAD_NODE_POINT\n');
    fprintf(foutload,'$ Generated by genPointLoads on %s\n', date);
    fprintf(foutload,'$ Normalization File:\n');
    fprintf(foutload,'$ %s\n', NormName);
    fprintf(foutload,'$ Normalization Isppa = %.1f W/cm^2\n', IsppaNorm);
end
