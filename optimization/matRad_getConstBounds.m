function [cl,cu] = matRad_getConstBounds(constraint,param,numOfScenarios,numOfVOIShiftScenarios)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT get constraint bounds function
% 
% call
%   [cl,cu] = matRad_getConstBounds(constraint,param)
%
% input
%   constraint: matRad constraint struct
%   param:      reference parameter
%
% output
%   cl: lower bounds on constraints
%   cu: lower bounds on constraints
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isequal(constraint.type, 'max dose constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'min dose constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'min mean dose constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'max mean dose constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'min max mean dose constraint') 

    cl = param(1);
    cu = param(2);

elseif isequal(constraint.type, 'min EUD constraint') 

    cl = param;
    cu = inf;

elseif isequal(constraint.type, 'max EUD constraint') 

    cl = -inf;
    cu = param;

elseif isequal(constraint.type, 'min max EUD constraint') 

    cl = param(1);
    cu = param(2);

elseif isequal(constraint.type, 'max DVH constraint') 

    cl = -inf;
    cu = constraint.volume/100;

    % alternative constraint calculation 1/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 1/4 %

elseif isequal(constraint.type, 'min DVH constraint') 

    cl = constraint.volume/100;
    cu = inf;

    % alternative constraint calculation 2/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 2/4 %
    
elseif isequal(constraint.type, 'max DCH constraint')
    
    cl = -inf;
    cu = constraint.coverage/100;
    
elseif isequal(constraint.type, 'min DCH constraint')
    
    cl = constraint.coverage/100;
    cu = inf;

elseif isequal(constraint.type, 'max DCH constraint2') || ...
       isequal(constraint.type, 'min DCH constraint2') 
   
    cl = -inf;
    cu = 0;
    
elseif isequal(constraint.type, 'max DCH constraint3')
    
    cl = -inf;
    cu = constraint.coverage/100;
    
elseif isequal(constraint.type, 'min DCH constraint3')  
    
    cl = constraint.coverage/100;
    cu = inf;
%     cl = [1 -inf -inf -inf -inf -inf -inf -inf -inf -inf]';
%     cu = inf(10,1);

elseif isequal(constraint.type, 'max DCH constraint4')
    
    cl = -inf(numOfScenarios,1);
    cu = repmat(constraint.volume/100,numOfScenarios,1);
    
elseif isequal(constraint.type, 'min DCH constraint4')  
    
    cl = repmat(constraint.volume/100,numOfScenarios,1);
    cu = inf(numOfScenarios,1);
    
elseif isequal(constraint.type, 'max DCH constraint5')
    
    cl = -inf(numOfVOIShiftScenarios,1);
    cu = repmat(constraint.volume/100,numOfVOIShiftScenarios,1);
    
elseif isequal(constraint.type, 'min DCH constraint5')  
    
    cl = repmat(constraint.volume/100,numOfVOIShiftScenarios,1);
    cu = inf(numOfVOIShiftScenarios,1);      

end % constraint switch
