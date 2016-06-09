function [dvhPoints,volume] = matRad_calcPDVH(coverage,doseVec,cst,numOfScenarios)

% set dose points in dvh
dvhPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);

% calculate DVH in every scenario
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
    doseInVoi   = doseVec{Scen}(cst{1,4}{1});
    dvh(Scen,:) = sum(bsxfun(@ge,doseInVoi,dvhPoints))/numel(cst{1,4}{1})*100;    
    end
elseif length(doseVec) == 1
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
        if isequal(cst{1,5}.VOIShift.shiftType,'rounded')
            dvh(Scen,:) = sum(bsxfun(@ge,doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.roundedShift.idxShift(Scen)),dvhPoints))/numel(cst{1,4}{1})*100;
            
        elseif isequal(cst{1,5}.VOIShift.shiftType,'linInterp')
            % lin interpolation in x
            c00 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y0Z0(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y0Z0(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c01 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y0Z1(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y0Z1(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c10 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y1Z0(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y1Z0(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c11 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y1Z1(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y1Z1(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
             
            % lin interpolation in y  
            c0  = c00.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen))+c10.*cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen);
            c1  = c01.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen))+c11.*cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen);
                
            % lin interpolation in z
            doseVecInterp = c0.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.z(Scen))+c1.*cst{1,5}.VOIShift.linInterpShift.idxShift.z(Scen);
            
            dvh(Scen,:) = sum(bsxfun(@ge,doseVecInterp,dvhPoints))/numel(cst{1,4}{1})*100;
        end

    end
end

% calculate PDVH
for j = 1:length(dvhPoints)
    VolumePointsSorted = sort(dvh(:,j),'descend');
    ix                 = max([1 ceil(coverage*numel(VolumePointsSorted))]);
    volume(1,j)        = VolumePointsSorted(ix);
end



end