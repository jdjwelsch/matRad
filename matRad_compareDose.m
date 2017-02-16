function matRad_compareDose(ct,stf,pln,cst,resultGUI,multScen)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compare external dose calculation to matRad dose calculation
% for same plan
% 
% call
%    
%
% input
%   cube1:          dose cube as an M x N x O array
%   cube2:          dose cube as an M x N x O array
%   resolution:     resolution of the cubes [mm/voxel]
%   criteria:      optional [1x2] vector depicting the gamma criteria                 
%
% output 
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% recalculate Dose
recalc_resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w,multScen);

% calculate gamma test
 dist2AgreeMm     = 3; % in [mm]
 relDoseThreshold = 3; % in [%]
 criteria = [dist2AgreeMm, relDoseThreshold]
[gC, mycm, gPassr] = matRad_gammaIndex(recalc_resultGUI.physicalDose, ...
    resultGUI.physicalDose, ct.resolution, criteria);

% TODO: plot dose slices along all axes for both distributions

% TODO: plot results gamma test

% (TODO: calculate and compare DVH for both distributions)

% TODO: adjust function to support biological effect comparison as well