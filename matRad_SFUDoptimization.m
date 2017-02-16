function [resultGUI] = matRad_SFUDoptimization(dij, cst, pln)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculation of single field uniform dose (SFUD) optimization

% call
%   [resultGUI] = SFUD_optimization(dij, cst, pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   (info:       struct containing information about optimization)
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust cst for single beams
sb_cst = cst;
for i=1:size(sb_cst,1)
    for j = 1:size(sb_cst{i,6},1)
        % biological dose splitting
        if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
            ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
            sb_cst{i,6}(j).dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                sb_cst{i,6}(j).dose/pln.numOfBeams +(sb_cst{i,6}(j).dose + ab));
        % physical dose splitting
        else
            sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
        end
    end
end


% initialise total weight vector
wTot = zeros(dij.totalNumOfBixels,1);

for i = 1:pln.numOfBeams
    % columns in total dij for single beam
    sb_col = find(dij.beamNum == i);
    % construct dij for single beam
    sb_dij.numOfBeams = 1;
    sb_dij.numOfVoxels = dij.numOfVoxels;
    sb_dij.resolution = dij.resolution;
    sb_dij.numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
    sb_dij.totalNumOfRays = sb_dij.numOfRaysPerBeam;
    sb_dij.totalNumOfBixels = size(sb_col, 1);
    sb_dij.dimensions = dij.dimensions;
    sb_dij.numOfScenarios = dij.numOfScenarios;
    sb_dij.ScenProb = dij.ScenProb;
    sb_dij.bixelNum = dij.bixelNum(sb_col);
    sb_dij.rayNum = dij.rayNum(sb_col);
    sb_dij.beamNum = dij.beamNum(sb_col);
    sb_dij.physicalDose{1} = dij.physicalDose{1}(:, sb_col);
    sb_dij.indexforOpt = dij.indexforOpt;
    if isfield(dij, 'RBE')
        sb_dij.RBE = dij.RBE;
    end
    if isfield(dij, 'mLETDose')
        sb_dij.mLETDose = dij.mLETDose(:, sb_col);
    end
    if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
        sb_dij.mAlphaDose{1} = dij.mAlphaDose{1}(:, sb_col);
        sb_dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:, sb_col);
    end
    
    % adjust pln to one beam only
    sb_pln = pln;
    sb_pln.gantryAngles = pln.gantryAngles(i);
    sb_pln.couchAngles = pln.couchAngles(i);
    
    % optimize single beam
    sb_resultGUI = matRad_fluenceOptimization(sb_dij,sb_cst,sb_pln);    
   
    % merge single beam weights into total weight vector
    wTot(sb_col) = sb_resultGUI.w;
end

% calculate dose
resultGUI = matRad_calcCubes(wTot,dij,cst,1);

end

