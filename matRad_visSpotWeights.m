function matRad_visSpotWeights(ct,stf,pln,cst,multScen,dij,resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to visualize Spot weights for single beams
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
%        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if  strcmp(pln.radiationMode, 'photons')
    fig = figure;
    % maximum weight
    w_max = max(resultGUI.w);
    for i=1:pln.numOfBeams
       rayPos_mat = vertcat(stf(i).ray(:).rayPos_bev);
       x_min = min(rayPos_mat(:, 1));
       z_min = min(rayPos_mat(:, 3));
       x_max = max(rayPos_mat(:, 1));
       z_max = max(rayPos_mat(:, 3));
       
       % initialise weight matrix
       weight_matrix = zeros((x_max-x_min)/pln.bixelWidth,(z_max-z_min)/pln.bixelWidth);
       
       for j=1:stf(i).numOfRays
           % find weight for this ray
           weightix = intersect(find(dij.beamNum==i), find(dij.rayNum==j));
           w_norm = sum(resultGUI.w(weightix))/length(weightix)/w_max; % normalize weight
           % find position in matrix corresponding to ray
           xpos = (stf(i).ray(j).rayPos_bev(1) + abs(x_min))/pln.bixelWidth +1;
           zpos = (stf(i).ray(j).rayPos_bev(3) + abs(z_min))/pln.bixelWidth +1;
           weight_matrix(xpos,zpos) = w_norm;          
       end
       
       % plot weight matrix
       ax = subplot(pln.numOfBeams,1, i);
       hold on;
       imagesc( [x_min x_max], [z_min z_max], weight_matrix);     
       title(sprintf(['spotweights in beam ' num2str(i)]));
       xlabel('x [mm]');
       ylabel('z [mm]');
       cmap = colormap(ax, 'jet');
       cmap(1,:) = [1 1 1];
       colormap(ax, cmap);
       cb = colorbar;
       ylabel(cb, 'normalized weight');
       hold off;
   end 
else
    % for protons and carbon Ions consider different energies as well
    
    % maximum weight
    w_max = max(resultGUI.w);
    for i=1:pln.numOfBeams
       rayPos_mat = vertcat(stf(i).ray(:).rayPos_bev);
       x_min = min(rayPos_mat(:, 1));
       z_min = min(rayPos_mat(:, 3));
       x_max = max(rayPos_mat(:, 1));
       z_max = max(rayPos_mat(:, 3));
       
       % find all energies used
       all_energies = unique(cat(2, stf(i).ray(:).energy));
       numOfEnergies = length(all_energies);
       % initialise weight matrix
       weight_matrix = zeros((x_max-x_min)/pln.bixelWidth,(z_max-z_min)/pln.bixelWidth, numOfEnergies);
       for j=1:numOfEnergies
           for k=1:stf(i).numOfRays
               weightix = intersect( intersect(find(dij.beamNum==i), find(dij.rayNum==k)), ... % find weights for this ray
                   find(stf(i).ray(k).energy(dij.bixelNum)==all_energies(j)));                 % find weight for this energy               
               
               w_norm = resultGUI.w(weightix)/w_max; % normalize weight
               % find position in matrix corresponding to ray
               xpos = (stf(i).ray(k).rayPos_bev(1) + abs(x_min))/pln.bixelWidth +1;
               zpos = (stf(i).ray(k).rayPos_bev(3) + abs(z_min))/pln.bixelWidth +1;
               weight_matrix(xpos,zpos, j) = w_norm;          
           end
           % plot weight matrix
%            plotsz1 = round(numOfEnergies/3); % split plots equally horizontally and vertically
%            plotsz2 = numOfEnergies- plotsz1;
%            ax = subplot(plotsz1,plotsz2, j);           
       end
       
       fig = figure('WindowScrollWheelFcn',@figScroll); % new figure for each beam
       ax = axes;
       curr_ix_energy = 1;
       img = imagesc( [x_min x_max], [z_min z_max], weight_matrix(:,:,curr_ix_energy));     
       title(sprintf(['spotweights in beam ' num2str(i) ', energy: ' num2str(all_energies(curr_ix_energy))]));
       xlabel('x [mm]');
       ylabel('z [mm]');
       cmap = colormap(ax, 'jet');
       cmap(1,:) = [1 1 1];
       colormap(ax, cmap);
       cb = colorbar;
       ylabel(cb, 'normalized weight');
       hold off;  

    end
    
end

%     function figScroll(src,callbackdata)
%       if (callbackdata.VerticalScrollCount > 0) && (curr_ix_energy<numOfEnergies)
%           curr_ix_energy = curr_ix_energy + 1;
%       elseif (callbackdata.VerticalScrollCount < 0) && (curr_ix_energy>1)
%           curr_ix_energy = curr_ix_energy - 1;
%       end
% 
%       img.CData = weight_matrix(:,:,curr_ix_energy);
%     end
