function [InputIntensity] = excludeFaceNodes(InputIntensity, mpn)
% function [InputIntensity] = excludeFaceNodes(InputIntensity, mpn)
%
% Exclude node intensity near the transducer face that may have violated the
% far-field assumption.
% 
% INPUTS:   InputIntensity - vector of intensity
%           mpn - node ID, x, y, z
%
% OUTPUTS:  InputIntensity - save as input, but with nodes near face zet to 0
%

AxialFaceTolCM = 0.001;

z0 = find(abs(mpn(:,4)) < AxialFaceTolCM);

InputIntensity(z0) = 0;
