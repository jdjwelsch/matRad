function matRad_visSpotWeights(stf,pln,dij,resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to visualize Spot weights for single beams
% 
% call
%    
%
% input
%   stf:              matRad stf struct
%   pln:              matRad pln struct
%   dij:              matRad dij struct
%   resultGUI:        matRad result struct for external calculation

% output 
%      plots visualisations for all beams - scroll to access different energies  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


counter = 0;
for i = 1:size(stf,2)
    for j = 1:stf(i).numOfRays
      for k = 1:stf(i).ray(j).numOfbIxelPerRay
          counter = counter + 1;
          
    end
end
    
    
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
           weight_matrix{i} = zeros((x_max-x_min)/pln.bixelWidth,(z_max-z_min)/pln.bixelWidth);

           for j=1:stf(i).numOfRays
               % find weight for this ray
               weightix = intersect(find(dij.beamNum==i), find(dij.rayNum==j));
               w_norm = sum(resultGUI.w(weightix))/length(weightix)/w_max; % normalize weight
               % find position in matrix corresponding to ray
               xpos = (stf(i).ray(j).rayPos_bev(1) + abs(x_min))/pln.bixelWidth +1;
               zpos = (stf(i).ray(j).rayPos_bev(3) + abs(z_min))/pln.bixelWidth +1;
               weight_matrix{i}(xpos,zpos) = w_norm;          
           end

           % plot weight matrix
           ax = subplot(pln.numOfBeams,1, i);
           hold on;
           imagesc( [x_min x_max], [z_min z_max], weight_matrix{i});     
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
           fprintf(['Calculate weights for Beam ' num2str(i) '...\n']);
            
           rayPos_mat = vertcat(stf(i).ray(:).rayPos_bev);
           x_min = min(rayPos_mat(:, 1));
           z_min = min(rayPos_mat(:, 3));
           x_max = max(rayPos_mat(:, 1));
           z_max = max(rayPos_mat(:, 3));

           % find all energies used
           all_energies = unique(cat(2, stf(i).ray(:).energy));
           numOfEnergies = length(all_energies);
           
           % initialise weight matrix
           weight_matrix{i} = zeros((x_max-x_min)/pln.bixelWidth,(z_max-z_min)/pln.bixelWidth, numOfEnergies);
           for j=1:numOfEnergies
               for k=1:stf(i).numOfRays
                   ray_ix = intersect(find(dij.beamNum==i), find(dij.rayNum==k));          % find weights for this ray
                   weightix = intersect(ray_ix, ...                                        % find weight for this energy               
                       min(ray_ix) + find(stf(i).ray(k).energy(dij.bixelNum(ray_ix))==all_energies(j)));

                   % normalize weight
                   w_norm = resultGUI.w(weightix)/w_max;
                   if isempty(w_norm)
                       w_norm = 0;
                   end
                   % find position in matrix corresponding to ray
                   xpos = (stf(i).ray(k).rayPos_bev(1) + abs(x_min))/pln.bixelWidth +1;
                   zpos = (stf(i).ray(k).rayPos_bev(3) + abs(z_min))/pln.bixelWidth +1;
                   weight_matrix{i}(xpos,zpos, j) = w_norm;          
               end
       
           end
           % new figure for each beam
           f = figure('Name', sprintf(['Beam ' num2str(i)]), 'WindowScrollWheelFcn', @figScroll);
           ax = axes;
           
           curr_ix_energy(i) = 1; % current energy slice
           
           % plot spot weights
           img = imagesc( [x_min x_max], [z_min z_max], weight_matrix{i}(:,:,curr_ix_energy(i)));     
           title(sprintf(['spotweights in beam ' num2str(i) ', energy: ' num2str(all_energies(curr_ix_energy(i)))]));
           xlabel('x [mm]');
           ylabel('z [mm]');
           cmap = colormap(ax, 'jet');
           %cmap(1,:) = [1 1 1]; % set color for 0 to white --> includes a
           %whole bin!
           caxis(ax, [0 1]); % set general range for colorbar
           colormap(ax, cmap);
           cb = colorbar;
           ylabel(cb, 'normalized weight');
           hold off;  

        end

    end


    function figScroll(src,callbackdata)
      % get handles for current figure and object
      curr_fig = gcf;
      h = gco(curr_fig);
      % get current beam number
      curr_beam = str2num(curr_fig.Name(end));
      
      % set new energy slice
      curr_ix_energy(curr_beam) = curr_ix_energy(curr_beam) - callbackdata.VerticalScrollCount;
      
      % project to allowed range
      curr_ix_energy(curr_beam) = min(curr_ix_energy(curr_beam), numOfEnergies);
      curr_ix_energy(curr_beam) = max(curr_ix_energy(curr_beam), 1);
      
      % change Data in plot
      h.CData = weight_matrix{curr_beam}(:,:,curr_ix_energy(curr_beam));
      title(sprintf(['spotweights in beam ' num2str(curr_beam) ', energy: ' num2str(all_energies(curr_ix_energy(curr_beam)))]));
      drawnow;
      
    end
end
