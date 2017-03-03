function sb_resultGUI = matRad_singleBeamView_nodij(ct,stf,pln,cst,multScen,weights,viewBeamNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  calculate single beam dose distributions without a dij struct being
%  present

% call
%   sb_resultGUI_nodij = matRad_singleBeamView_nodij(pln, dij, cst, resultGUI, viewBeamNum)
%
% input
%   cst:         matRad cst struct
%   weights:     weights for beamspots (usually matRad resultGUI.w struct)
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
    sb_resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,weights,multScen);

elseif ~isnumeric(viewBeamNum)
    fprintf('Error: wrong beam number format in singleBeamView \n');

else
    % calculate single beam dose
    
    % adjust cst for single beams
    sb_cst = cst;
    for i=1:size(sb_cst,1)
        for j = 1:size(sb_cst{i,6},1)
            % biological dose splitting
            if strcmp(pln.bioOptimization, 'LEMIV_effect') || strcmp(pln.bioOptimization, 'LEMIV_RBExD')
                ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
                sb_cst{i,6}(j).dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                    sb_cst{i,6}(j).dose/pln.numOfBeams +(sb_cst{i,6}(j).dose + ab));
            % physical dose splitting
            else
                sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
            end
        end
    end


    fprintf(['Calculating Dose Beam ' num2str(viewBeamNum) '\n']);

    % singlebeam stf
    sb_stf = stf(viewBeamNum);

    % singlebeam weights
    if viewBeamNum==1
        offset = 0;
    else
        offset = sum(stf(1:(viewBeamNum-1)).totalNumOfBixels);
    end
    sb_weights = weights(offset + (1:stf(viewBeamNum).totalNumOfBixels));
    
    % adjust pln to one beam only
    sb_pln = pln;
    sb_pln.isoCenter = pln.isoCenter(viewBeamNum,:);
    sb_pln.numOfBeams = 1;
    sb_pln.gantryAngles = pln.gantryAngles(viewBeamNum);
    sb_pln.couchAngles = pln.couchAngles(viewBeamNum);

    % calculate single beam dose
    sb_resultGUI = matRad_calcDoseDirect(ct,sb_stf,sb_pln,sb_cst,sb_weights,multScen);

    % keep full set of weights and for other beams
    % sb_resultGUI.w = weights;

end

end