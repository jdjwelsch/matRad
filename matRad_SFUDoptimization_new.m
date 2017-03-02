function [resultGUI] = matRad_SFUDoptimization_new(ct,stf,pln,cst,multScen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   [1]    https://ro-journal.biomedcentral.com/articles/10.1186/s13014-016-0705-8
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust cst for single beams
sb_cst = cst;
for i=1:size(sb_cst,1)
    for j = 1:size(sb_cst{i,6},1)
        % biological dose splitting
        if strcmp(pln.bioOptimization, 'LEMIV_effect') || strcmp(pln.bioOptimization, 'LEMIV_RBExD')
            ab = sb_cst{i,5}.alphaX / sb_cst{i,5}.betaX;
            % dose per fraction
            fx_dose = sb_cst{i,6}(j).dose / pln.numOfFractions;
            % calculate dose per beam per fraction according to [1]
            fx_dose = -0.5*ab +sqrt( 0.25*ab^2 + ...
                fx_dose/pln.numOfBeams *(fx_dose + ab));
            % calculate pseudo total Dose per Beam
            sb_cst{i,6}(j).dose = fx_dose * pln.numOfFractions;
            
        % physical dose splitting
        else
            sb_cst{i,6}(j).dose = sb_cst{i,6}(j).dose/pln.numOfBeams;
        end
    end
end

% initialise total weight vector
wTot = [];

for i = 1:pln.numOfBeams
    fprintf(['optimizing beam ' num2str(i) '...\n']);
    % single beam stf
    sb_stf = stf(i);
    
    % adjust pln to one beam only
    sb_pln = pln;
    sb_pln.isoCenter = pln.isoCenter(i,:);
    sb_pln.numOfBeams = 1;
    sb_pln.gantryAngles = pln.gantryAngles(i);
    sb_pln.couchAngles = pln.couchAngles(i);
    
    % calculate single beam dij
    sb_dij = matRad_calcParticleDose(ct,sb_stf,sb_pln,sb_cst,multScen,false);
    
    % optimize single beam
    sb_resultGUI = matRad_fluenceOptimization(sb_dij,sb_cst,sb_pln);    
   
    % merge single beam weights into total weight vector
    wTot = [wTot ; sb_resultGUI.w];
    
end

fprintf('Calculate total dose...\n');
% calculate dose
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,wTot,multScen);

end

