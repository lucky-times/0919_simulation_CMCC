function LTE_plot_capesso_files(eNodeBs,pathlossmaps, networkMacroscopicPathlossMap, capesso_params)
% Plots the loaded Capesso files based on the data stored on the struct
%
% (c) Josep Colom Ikuno, Martin Taranetz INTHFT, 2010

UEs = [];
map_res = capesso_params.maps_resolution;
current_TTI = 1;
global LTE_config;
first_figure = LTE_config.plots.capesso_maps_begin+0;
eNodeBs_pos = LTE_plot_show_network(eNodeBs, UEs, map_res, current_TTI,first_figure,'Loaded Capesso maps and ROI');
figure(first_figure);
hold on;
% Set a common pathloss scale and plot the pathloss files in the overall grid
for i_=1:length(pathlossmaps)
    if i_==1
        c_axis  = [min(pathlossmaps(i_).pathloss_map(:))  max(pathlossmaps(i_).pathloss_map(:))];
        c_axis2 = [min(pathlossmaps(i_).elevation_map(:)) max(pathlossmaps(i_).elevation_map(:))];
    else
        c_axis  = [min([c_axis(1);  pathlossmaps(i_).pathloss_map(:)])  max([c_axis(2);  pathlossmaps(i_).pathloss_map(:)])];
        c_axis2 = [min([c_axis2(1); pathlossmaps(i_).elevation_map(:)]) max([c_axis2(2); pathlossmaps(i_).elevation_map(:)])];
    end
    SWxmap = pathlossmaps(i_).description.SWxmap; % Lower-leftmost corner (SW)
    SWymap = pathlossmaps(i_).description.SWymap; % Lower-leftmost corner (SW)
    height = pathlossmaps(i_).description.ncols * pathlossmaps(i_).description.ydim;
    width  = pathlossmaps(i_).description.nrows * pathlossmaps(i_).description.xdim;
    rectangle('Position',[SWxmap, SWymap, width, height]);
    text(SWxmap,SWymap,num2str(i_),'FontSize',10);
end
xlimits = xlim;
ylimits = ylim;
ROI_increase_factor = 0.05;
xlimits = ROI_increase_factor*abs(xlimits(1)-xlimits(2))*[-1 1] + xlimits;
ylimits = ROI_increase_factor*abs(ylimits(1)-ylimits(2))*[-1 1] + ylimits;
xlim(xlimits);
ylim(ylimits);

% Draw a suggested ROI
sug_roi_min    = min(eNodeBs_pos);
sug_roi_max    = max(eNodeBs_pos);
sug_roi_width  = sug_roi_max(1) - sug_roi_min(1);
sug_roi_height = sug_roi_max(2) - sug_roi_min(2);
eNodeB_ROI_increase_factor = capesso_params.eNodeB_ROI_increase_factor;
rectangle('Position',[sug_roi_min(1)-eNodeB_ROI_increase_factor*sug_roi_width, sug_roi_min(2)-eNodeB_ROI_increase_factor*sug_roi_height, sug_roi_width*(1+2*eNodeB_ROI_increase_factor), sug_roi_height*(1+2*eNodeB_ROI_increase_factor)],'EdgeColor','red');

hold off;

% Plot pathloss for each capesso file
pathloss_figure_idx = LTE_config.plots.capesso_maps_begin+1;
for i_=1:length(pathlossmaps)
    pathlossmap = pathlossmaps(i_);
    cleaned_filename = strrep(pathlossmap.description.filename,'#2F','/');
    
    % Plot pathloss
    figure_steps_col = 5;
    figure_steps_row = 5;
    pncols_ = pathlossmap.description.ncols;
    pnrows_ = pathlossmap.description.nrows;
    pxdim_  = pathlossmap.description.xdim;
    pydim_  = pathlossmap.description.ydim;
    pathloss_colormap = [1 0 0;1 0.06275 0;1 0.1294 0;1 0.1922 0;1 0.2588 0;1 0.3216 0;1 0.3882 0;1 0.451 0;1 0.5176 0;1 0.5804 0;1 0.6471 0;1 0.7098 0;1 0.7725 0;1 0.8392 0;1 0.902 0;1 0.9686 0;0.9686 1 0;0.902 1 0;0.8392 1 0;0.7725 1 0;0.7098 1 0;0.6471 1 0;0.5804 1 0;0.5176 1 0;0.451 1 0;0.3882 1 0;0.3216 1 0;0.2588 1 0;0.1922 1 0;0.1294 1 0;0.06275 1 0;0 1 0;0 1 0.06275;0 1 0.1255;0 1 0.1882;0 1 0.251;0 1 0.3137;0 1 0.3765;0 1 0.4392;0 1 0.502;0 1 0.5608;0 1 0.6235;0 1 0.6863;0 1 0.749;0 1 0.8118;0 1 0.8745;0 1 0.9373;0 1 1;0 0.9373 1;0 0.8745 1;0 0.8118 1;0 0.749 1;0 0.6863 1;0 0.6235 1;0 0.5608 1;0 0.498 1;0 0.4392 1;0 0.3765 1;0 0.3137 1;0 0.251 1;0 0.1882 1;0 0.1255 1;0 0.06275 1;0 0 1];
    xTicks  = 0 : pncols_/figure_steps_col: pncols_;
    xLabels = 0:pncols_*pxdim_/figure_steps_col:pncols_*pxdim_;
    yTicks  = 0 : pnrows_/figure_steps_row: pnrows_;
    yLabels = 0:pnrows_*pydim_/figure_steps_row:pnrows_*pydim_;
    
    figure_capesso_los_file = figure(pathloss_figure_idx+i_-1);
    colormap('default');
    
    % Pathloss map (1st subplot)
    axes1 = subplot(3+capesso_params.plot_antenna_gain_patterns ,3,1,'Parent',figure_capesso_los_file);
    x_axis = [0 pncols_*pxdim_] + pathlossmap.description.NWxmap;
    y_axis = [0 pnrows_*pydim_] + pathlossmap.description.NWymap;
    %imagesc(x_axis,y_axis,-pathlossmap.pathloss_map);
    imagesc(pathlossmap.description.roi_x,pathlossmap.description.roi_y,-pathlossmap.pathloss_map);
    set(gca,'YDir','normal');
    caxis(axes1,-[c_axis(2) c_axis(1)]);
    colorbar;
    grid(axes1,'on');
    title(axes1,sprintf('%d-Capesso: %s ',i_,strrep(cleaned_filename,'_','\_')));
    hold(axes1,'on');
    for b_=1:length(eNodeBs)
        [faces verts] = get_patch_data(eNodeBs(b_),pathlossmap.description.xdim);
        patch('Faces',faces,'Vertices',verts,'EdgeColor','k','LineWidth',2);
        %set(p,'EdgeColor','none');
    end
    xlim(axes1,networkMacroscopicPathlossMap.roi_x);
    ylim(axes1,networkMacroscopicPathlossMap.roi_y);
    hold(axes1,'off');
    
    % Elevation map (2nd subplot)
    axes2 = subplot(3+capesso_params.plot_antenna_gain_patterns ,3,2,'Parent',figure_capesso_los_file);
    imagesc(pathlossmap.description.roi_x,pathlossmap.description.roi_y,pathlossmap.elevation_map);
    set(gca,'YDir','normal');
    caxis(axes2,[c_axis2(1) c_axis2(2)]);
    grid(axes2,'on');
    hold(axes2,'on');
    for b_=1:length(eNodeBs)
        [faces verts] = get_patch_data(eNodeBs(b_),pathlossmap.description.xdim);
        patch('Faces',faces,'Vertices',verts,'EdgeColor','k','LineWidth',2);
        %set(p,'EdgeColor','none');
    end
    hold(axes2,'off');
    xlim(axes2,networkMacroscopicPathlossMap.roi_x);
    ylim(axes2,networkMacroscopicPathlossMap.roi_y);
    xlim(axes2,networkMacroscopicPathlossMap.roi_x);
    ylim(axes2,networkMacroscopicPathlossMap.roi_y);
    colorbar;
    title(axes2,sprintf('%d-Elevation map (m). %s area',i_,strrep(cleaned_filename,'_','\_')));
    
    % eNodeB position (3rd subplot)
    axes3 = subplot(3+capesso_params.plot_antenna_gain_patterns,3,3,'Parent',figure_capesso_los_file);
    hold(axes3,'on');
    for b_=1:length(eNodeBs)
        % Plot a line that tells where the antennas are pointing
        vector_length = 80;
        origin = eNodeBs(b_).pos;
        for s_=1:length(eNodeBs(b_).sectors)
            angle = wrapTo360(-eNodeBs(b_).sectors(s_).azimuth+90);
            vector = vector_length*[ cosd(angle) sind(angle) ];
            destiny = vector + origin;
            
            plot([origin(1) destiny(1)],[origin(2) destiny(2)]);
        end
        % Plot the eNodeBs
        scatter(axes3,eNodeBs(b_).pos(1),eNodeBs(b_).pos(2),'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','black');
        text(eNodeBs(b_).pos(1)+map_res*5,eNodeBs(b_).pos(2),num2str(eNodeBs(b_).id));
    end
    hold(axes3,'off');
    grid(axes3,'on');
    xlim(axes3,networkMacroscopicPathlossMap.roi_x);
    ylim(axes3,networkMacroscopicPathlossMap.roi_y);
    colorbar; % This is just so the size of the plot is the same as the one above
    title(axes3,sprintf('eNodeB posititions'));
    
    % Plot sector antenna gain (assuming 3 sectors!!)
    % and pathlossmap with applied antenna gain in ROI for every sector

    for s_ = 1:length(eNodeBs(b_).sectors)
        
        % Plot antenna gain maps
        azimuth = eNodeBs(i_).sectors(s_).azimuth;
        current_axes = subplot(3+capesso_params.plot_antenna_gain_patterns,3,3+s_);
        imagesc(pathlossmap.description.roi_x,pathlossmap.description.roi_y,pathlossmap.antenna_gain_map(:,:,s_));
        set(gca,'YDir','normal');
        hold(current_axes,'on');
        for b_=1:length(eNodeBs)
            [faces verts] = get_patch_data(eNodeBs(b_),pathlossmap.description.xdim);
            patch('Faces',faces,'Vertices',verts,'EdgeColor','k','LineWidth',2);
            %set(p,'EdgeColor','none');
        end
        hold(current_axes,'off');
        colorbar;
        title(strrep(sprintf('Antenna Gain (Sec%d), \\theta=%3.0f°',s_, azimuth),'_','\_'));
        xlim(networkMacroscopicPathlossMap.roi_x);
        ylim(networkMacroscopicPathlossMap.roi_y);
        caxis([-60,20]);
        
        % Plot pathlossmaps with applied antenna gain in ROI
        axes1 = subplot(3+capesso_params.plot_antenna_gain_patterns,3,6+s_,'Parent',figure_capesso_los_file);
        imagesc(networkMacroscopicPathlossMap.roi_x, networkMacroscopicPathlossMap.roi_y, -networkMacroscopicPathlossMap.pathloss(:,:,s_,i_));
        set(gca,'YDir','normal');
        caxis(axes1,-[c_axis(2) c_axis(1)]);
        colorbar;
        title(strrep(sprintf('Overlay in ROI (Sec%d)',s_),'_','\_'));
        
        % Plot antenna gain patterns using polar plot for each sector without
        % respecting the mechanical downtilt
        if capesso_params.plot_antenna_gain_patterns > 0
            sec_electrical_downtilt = eNodeBs(i_).sectors(s_).electrical_downtilt;
            sec_antenna             = eNodeBs(i_).sectors(s_).antenna;
            
            if find(sec_antenna.electrical_tilt == sec_electrical_downtilt)
                index_ = find(sec_antenna.electrical_tilt == sec_electrical_downtilt, 1, 'first');
            else
                error('Gain pattern for electrical tilt of %f° not available !\n', sec_electrical_downtilt);
                index_ = 0;
            end
            
            if index_  > 0
                %             axes4 = subplot(3+capesso_params.plot_antenna_gain_patterns,3,9+s_,'Parent',figure_capesso_los_file);
                %             polar2((0:359)'/180*pi,  10.^((sec_antenna.max_antenna_gain(index_)-sec_antenna.horizontal_gain_pattern(:,index_))/10)); hold on;
                %             polar2(-(0:359)'/180*pi, 10.^((sec_antenna.max_antenna_gain(index_)-sec_antenna.vertical_gain_pattern(:,index_))/10), 'red');
                %             legend('horizontal [dB]', 'vertical [dB]', 2,'Location','SouthEastOutside');
                %             title(['Antenna Type ' eNodeBs(i_).sectors(s_).antenna_type ' (Electrical Tilt ',num2str(sec_antenna.electrical_tilt(index_)),'°)']);
                %             hold off;
                
                axes4 = subplot(3+capesso_params.plot_antenna_gain_patterns,3,9+s_,'Parent',figure_capesso_los_file);
                plot_tilt = sec_electrical_downtilt; % Plot with the electrical tilt from the sector. The effect of the mechanical tilt is not shown in the polar plot, but it is applied when calculating the pathloss map.
                data_limits = [-15 0 3];
                [hor_degrees hor_gain ver_degrees ver_gain max_gain] = sec_antenna.gain_patterns(plot_tilt);
                utils.miscUtils.polar2(hor_degrees/180*pi, hor_gain-max_gain, data_limits,'blue');
                hold(axes4,'all');
                utils.miscUtils.polar2(ver_degrees/180*pi, ver_gain-max_gain, data_limits,'red');
                title(axes4,sprintf('%d antenna. %3.0f° electrical tilt',sec_antenna.antenna_type,plot_tilt));
                hold(axes4,'off');
            end
            
        end
    end
end

function [faces verts] = get_patch_data(an_eNodeB,data_res)
% Processes the position data from an eNodeB to generate data which can be fed to the "patch" function.
num_sectors = length(an_eNodeB.sectors);
verts = zeros(num_sectors,2);
faces = zeros(num_sectors,2);
pos = an_eNodeB.pos;
vector_length = 1*data_res;

for s_=1:num_sectors
    angle = wrapTo360(-an_eNodeB.sectors(s_).azimuth+90);
    vector = vector_length*[ cosd(angle) sind(angle) ];
    verts(s_,:) = vector + pos;
    if s_>1
        faces(s_-1,:) = [s_-1 s_];
    end
end
faces(s_,:) = [s_ 1];
