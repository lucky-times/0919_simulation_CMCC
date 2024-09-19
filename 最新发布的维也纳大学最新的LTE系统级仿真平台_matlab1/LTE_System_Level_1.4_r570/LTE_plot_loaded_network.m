function LTE_plot_loaded_network(eNodeBs,networkPathlossMap,varargin)
% Function that shows a few plots after a network has been loaded from a file
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at

global LTE_config;

if ~isempty(varargin)
    plot_shadow_fading = true;
    networkShadowFadingMap = varargin{1};
else
    plot_shadow_fading = false;
end

% Not very clean, but will do the trick
roi_to_map_x = networkPathlossMap.roi_x;
roi_to_map_y = networkPathlossMap.roi_y;

%% Plot the antenna gain pattern
if LTE_config.show_network>0
    all_sectors = [eNodeBs.sectors];
    all_antennas = [all_sectors.antenna];
    [anntenna_types, m, n] = unique([all_antennas.antenna_type]);
    number_of_subplots_row_col = ceil(sqrt(length(anntenna_types)));
    
    % Plot antenna gain patter for each different antenna in the simulation
    figure(LTE_config.plots.antenna_gain_pattern);
    for ant_idx = 1:length(anntenna_types)
        example_ant_idx = find(n==ant_idx,1,'first');
        an_antenna = all_antennas(example_ant_idx);
        
        % Choose correct antenna pattern
        the_axes = subplot(number_of_subplots_row_col,number_of_subplots_row_col,ant_idx,'replace');
        
        if ~an_antenna.pattern_is_3D
            % 2D antenna
            angle = -180:0.1:180;
            gain = zeros(1,length(angle));
            for i_=1:length(angle)
                gain(i_) = an_antenna.gain(angle(i_));
            end
            cla(the_axes);
            plot(the_axes,angle,gain);
            ylim(the_axes,ylim*1.1);
            title(the_axes,{['Antenna gain, ' an_antenna.antenna_type ' antenna']});
            xlabel(the_axes,{'\theta [°]'});
            ylabel(the_axes,{'gain [dB]'});
            box(the_axes,'on');
            grid(the_axes,'on');
        else
            % 3D antenna
            plot_tilt = 0; % Plot with 0° electrical (and mechanical) tilt
            data_limits = [-15 0 3];
            [hor_degrees hor_gain ver_degrees ver_gain max_gain] = an_antenna.gain_patterns(plot_tilt);
            utils.miscUtils.polar2(hor_degrees/180*pi, hor_gain-max_gain, data_limits,'blue');
            hold(the_axes,'all');
            utils.miscUtils.polar2(ver_degrees/180*pi, ver_gain-max_gain, data_limits,'red');
            title(the_axes,sprintf('%d antenna.\nblue: hor, red: vert [dBi],\n%3.0f° electrical tilt',an_antenna.antenna_type,plot_tilt));
            hold(the_axes,'off');
        end
    end
end

%% Plot of how the macroscopic pathloss looks like
if LTE_config.show_network>0 && LTE_config.macroscopic_pathloss_is_model
    macroscopic_pathloss_model = LTE_common_get_macroscopic_pathloss_model;
    % Will set the maximum distance as the diagonal that crosses the ROI
    range = sqrt((roi_to_map_x(2)-roi_to_map_x(1))^2+(roi_to_map_y(2)-roi_to_map_y(1))^2);
    distances = 0:LTE_config.map_resolution:range;
    pathlosses = macroscopic_pathloss_model.pathloss(distances);
    figure(LTE_config.plots.macroscopic_pathloss);
    clf;
    plot(distances,pathlosses);
    title(['Macroscopic pathloss, using ' macroscopic_pathloss_model.name ' model']);
    xlabel('Distance [m]');
    ylabel('Pathloss [dB]');
    box on;
    grid on;
end

%% Plot of sector macroscopic pathlosses (ASSUMING 3 SECTORS)
if LTE_config.show_network>0
    for sector_idx = 1:3
        figure(LTE_config.plots.macroscopic_pathloss_sector(sector_idx));
        number_cols = 3;
        number_rows = ceil(length(eNodeBs)/number_cols);
        for b_ = 1:length(eNodeBs)
            subplot(number_rows,number_cols,b_);
            imagesc(networkPathlossMap.roi_x,networkPathlossMap.roi_y,networkPathlossMap.pathloss(:,:,sector_idx,b_));
            set(gca,'YDir','normal');
            title(sprintf('eNodeB %d sector %d',b_,sector_idx));
            colorbar;
            hold on;
            scatter(eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'MarkerEdgeColor','white','MarkerFaceColor','black');
            text(eNodeBs(b_).pos(1)+5*networkPathlossMap.data_res,eNodeBs(b_).pos(2),num2str(b_),'Color','w');
        end
    end
end

%% Plot shadow fading
if LTE_config.show_network>0 && plot_shadow_fading
    num_eNodeBs = length(eNodeBs);
    N_cols = 3;
    N_rows = ceil(num_eNodeBs/N_cols);
    figure(LTE_config.plots.shadow_fading_loss);
    for i_=1:num_eNodeBs
        subplot(N_rows,N_cols,i_);
        imagesc(networkShadowFadingMap.roi_x,networkShadowFadingMap.roi_y,networkShadowFadingMap.pathloss(:,:,i_));
        set(gca,'YDir','normal');
        title(['Shadow fading, eNodeB ' num2str(i_)]);
        colorbar;
    end
    
%     % Histogram, to see if they are really gaussian or not
%     size_xy = size(networkShadowFadingMap.pathloss,1)*size(networkShadowFadingMap.pathloss,2);
%     figure(LTE_config.plots.shadow_fading_loss_histogram);
%     for i_=1:num_eNodeBs
%         subplot(N_rows,N_cols,i_);
%         reshaped_map = reshape(networkShadowFadingMap.pathloss(:,:,i_),1,size_xy);
%         map_mean = mean(reshaped_map);
%         map_sd = std(reshaped_map);
%         hist(reshaped_map,75);
%         if i_==1
%             xlimits = xlim;
%             ylimits = ylim*1.2;
%         end
%         xlim(xlimits);
%         ylim(ylimits);
%         grid on;
%         set(gca,'YDir','normal');
%         title(['Shadow fading, eNodeB ' num2str(i_) ' mean: ' num2str(map_mean,'%3.2f') ' sd: ' num2str(map_sd,'%3.2f') ]);
%     end
end
