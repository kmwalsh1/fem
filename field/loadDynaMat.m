function [intensity, focal_depth, mpn] = loadDynaMat(infile)
% function [intensity, focal_depth, mpn] = loadDynaMat(infile)
% 
% Load dyna*.mat file that includes the intensity data and the focal
% configuration struct.
%
% INPUTS:   infile (str) - input file
%
% OUTPUTS:  intensity (float vector) - intensity data vector
%           focal_depth (float) - focal depth (cm)
%           mpn (# nodex x 4) - node ID, x, y, z
%

load(infile);

mpn = FIELD_PARAMS.mpn;

focal_depth = FIELD_PARAMS.focus(3)*100; % m -> cm
