function sb_resultGUI = matRad_singleBeamView(pln, dij, cst, resultGUI, viewBeamNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot single beam dose distributions

% call
%   resultGUI = matRad_singleBeamView(pln, dij, cst, resultGUI, viewBeamNum)
%
% input
%   dij:         matRad dij struct
%   cst:         matRad cst struct
%   resultGUI:   matRad resultGUI struct
%   viewBeamNum: beam to be visualised, either integer for beamNum or 'all'
%                to calculate total dose
%
% output
%   
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(viewBeamNum, 'all')
    % recalculate standard resultGUI for all beams
    fprintf('Calculating total dose \n');
    sb_resultGUI = matRad_calcCubes(resultGUI.w,dij,cst,1);

elseif ~isnumeric(viewBeamNum)
    fprintf('Error: wrong beam number format in sbView \n');

else
    % calculate single beam dose
    
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

    for i = viewBeamNum
        fprintf(['Calculating Dose Beam ' num2str(viewBeamNum) '\n']);
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

        sb_w = resultGUI.w(sb_col);
        sb_resultGUI = matRad_calcCubes(sb_w,sb_dij,sb_cst,1);
        
        % keep full set of weights and for other beams
        sb_resultGUI.w = resultGUI.w;
    end
end
end