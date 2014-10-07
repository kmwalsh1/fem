function checkLat0(mpn)
% function checkLat0(mpn)
%
% Check that the mesh has a lat = 0 position to make sure that the center of
% the ROE is properly sampled.
%
% INPUTS:   mpn - node ID, x, y, z matrix
%
% OUTPUTS:  Errors if lat = 0 does not exist.
%

LatTolCM = 1e3; % search tol (cm)

if(isempty(find(abs(mpn(:,3)) < LatTolCM))),
    error('lat = 0 nodes missing for rad force excitation');
end;
