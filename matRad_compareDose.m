function matRad_compareDose(ct,stf,pln,cst,multScen,resultGUI,recalc_resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compare external dose calculation to matRad dose calculation
% for the same plan
% 
% call
%    
%
% input
%   ct:               matRad ct struct
%   stf:              matRad stf struct
%   pln:              matRad pln struct
%   cst:              matRad cst cell
%   multScen:         matRad multScen struct
%   resultGUI:        matRad result struct for external calculation
%   recalc_resultGUI: optional matRad result struct - will be calculated if
%                     not provided

% output 
%    plots with slices through result cube, results of gamma Test,
%    dhvs for both distributions    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin<7
    % recalculate Dose if recalculation not provided
    recalc_resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w,multScen);
end
%% set start parameters
% consider biological effect comparison
if strcmp(pln.bioOptimization, 'RBExDose')
    result = resultGUI.RBExDose;
    recalc_result = recalc_resultGUI.RBExDose;
    optimizationQuantity = 'RBExDose [Gy]';
elseif strcmp(pln.bioOptimization, 'LEMIV_effect')
    result = resultGUI.effect;
    recalc_result = recalc_resultGUI.effect;
    optimizationQuantity = 'effect';
elseif strcmp(pln.bioOptimization, 'LEMIV_RBExDose')
    result = resultGUI.RBExDose;
    recalc_result = recalc_resultGUI.RBExDose;
    optimizationQuantity = 'RBExDose';
else    
    result = resultGUI.physicalDose;
    recalc_result = recalc_resultGUI.physicalDose;
    optimizationQuantity = 'Dose [Gy]';
end

% calculate gamma test
dist2AgreeMm     = 3; % in [mm]
relDoseThreshold = 3; % in [%]
resolution(1) = ct.resolution.x;
resolution(2) = ct.resolution.y;
resolution(3) = ct.resolution.z;

[gammaCube, myColormap, gammaPassRate] = matRad_gammaIndex(recalc_result, ...
    result, resolution, [dist2AgreeMm, relDoseThreshold]);


%% plot results

% dose slices along all axes through isocenter

% arrays for plotting
xarray = (1:ct.cubeDim(1))*resolution(1);
yarray = (1:ct.cubeDim(1))*resolution(2);
zarray = (1:ct.cubeDim(3))*resolution(3);
% isocenter indices in cubes
xiso = round(stf(1).isoCenter(1)/resolution(1));
yiso = round(stf(1).isoCenter(2)/resolution(2));
ziso = round(stf(1).isoCenter(3)/resolution(3));

figure;
subplot(2, 3, 1);
hold on
plot(xarray, result(:, yiso, ziso), ':r', 'LineWidth', 1.5);
plot(xarray, recalc_result(:, yiso, ziso), 'b', 'LineWidth', 1.5);
xlabel('x [mm]');
ylabel(optimizationQuantity);
title('slice through yiso, ziso');
legend('Syngo Dose', 'matRad Dose')
hold off

subplot(2, 3, 2);
hold on
plot(yarray, result(xiso, :, ziso), ':r', 'LineWidth', 1.5);
plot(yarray, recalc_result(xiso, :, ziso), 'b', 'LineWidth', 1.5);
xlabel('y [mm]');
ylabel(optimizationQuantity);
title('slice through xiso, ziso');
legend('Syngo Dose', 'matRad Dose')
hold off

subplot(2, 3, 3);
hold on
plot(zarray, squeeze(result(xiso, yiso, :)), ':r', 'LineWidth', 1.5);
plot(zarray, squeeze(recalc_result(xiso, yiso, :)), 'b', 'LineWidth', 1.5);
xlabel('z [mm]');
ylabel(optimizationQuantity);
title('slice through xiso, yiso');
legend('Syngo Dose', 'matRad Dose')
hold off

% plot difference in z=ziso plane
resdiff = abs(result-recalc_result);
subplot(2, 3, 4)
hold on
title('absolute difference in xy-plane at z=ziso');
imagesc(resdiff(:, :, ziso));
xlabel('x [mm]');
ylabel('y [mm]');
colormap('jet')
cb = colorbar;
ylabel(cb, optimizationQuantity);
hold off;

% plot results gamma test
subplot(2, 3, 5)
hold on;
imagesc(gammaCube(:,:,ziso),[0 2])
xlabel('x [mm]');
ylabel('y [mm]');
colormap(myColormap);
colorbar
title(['gamma pass rate ' num2str(gammaPassRate,5)  ' %, gamma criterion (' ...
    num2str(relDoseThreshold) '% / ' num2str(dist2AgreeMm) 'mm)']);
hold off;

% calculate and compare DVH for both distributions
figure;
fprintf('Calculating DVH for Syngo dose: \n');
matRad_calcDVH(resultGUI,cst,pln,1)
fprintf('Calculating DVH for matRad dose: \n');
matRad_calcDVH(recalc_resultGUI,cst,pln,2)

