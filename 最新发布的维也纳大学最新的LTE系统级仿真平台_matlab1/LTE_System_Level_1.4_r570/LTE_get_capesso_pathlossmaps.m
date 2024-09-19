function [pathlossmaps unused_pathloss_maps] = LTE_get_capesso_pathlossmaps(capesso_params,eNodeBs)
% Reads out data files from capesso. Outputs a matrix containing the
% pathloss maps (eNodeB_idx,sector_idx) in the same order as in the eNodeBs
% object matrix. For capesso file the matrix is of size
% (length(eNodeBs),1), as only one pathloss file per eNodeB is needed. The
% second output logs interpreted pathloss files that were not output (the
% eNodeB/site was not on the list or only one sector was needed, thus the
% others were discarded).
%
% (c) Martin Taranetz, Josep Colom Ikuno INTHFT, 2010

%% Parameters

folder = capesso_params.pathloss_data_folder;

% Where to get the digital terrain model
dtm_folder          = capesso_params.dtm_folder;
dtm_file_name       = capesso_params.dtm_file_name;
dtm_hdr_file_name   = capesso_params.dtm_hdr_file_name;
%enable_dtm_plotting = capesso_params.enable_debug_plotting;

% Enable / disable plots of the pathloss maps
enable_debug_plotting     = capesso_params.enable_debug_plotting;


%% Get full DTM map and cut out region of interest for every pathlossmap
% (here, the region of interest denotes the size of the Capesso
% pathlossmap)
elevation_map_full = LTE_init_generate_dtm(dtm_folder,dtm_file_name,dtm_hdr_file_name,false);

%% Read Capesso Files and generate pathloss map

% .clos files are renamed .zip files containing the .los files
% If there are any .clos files in the directory, they are shifted to the
% zip_files directory and renamed to .zip
% In a second step, the .zip files are extracted back in the folder of
% origin, which also contains the .par description files

clos_files = dir(fullfile(folder,'*.clos'));
if ~isempty(clos_files)
    mkdir(folder, 'zip_files');
    for i_=1:1:length(clos_files)
        zipfilename = [strtok(clos_files(i_).name,'.'),'.zip'];
        movefile(fullfile(folder,clos_files(i_).name),[folder,'/zip_files/',zipfilename]);
        unzip([folder,'/zip_files/',zipfilename], folder);
    end
end

% unzipped .los files are ready for readout
los_files = dir(fullfile(folder,'*.los'));

unused_pathloss_maps = [];

for i_=1:1:length(los_files)
    
    file_name = los_files(i_).name;
    plfile_ = fopen(fullfile(folder, file_name),'r');
    %Search for .par file and get description for .los file
    pathlossmap.description = get_pathloss_files_description(file_name, folder);
    
    PATHLOSS_MAP_A = fread(plfile_,'int16');
    PATHLOSS_MAP_A = reshape(PATHLOSS_MAP_A,pathlossmap.description.ncols, pathlossmap.description.nrows);

    PATHLOSS_MAP_B = PATHLOSS_MAP_A.';
    PATHLOSS_MAP_B = PATHLOSS_MAP_B./16; 

    pathlossmap.pathloss_map(:,:)= PATHLOSS_MAP_B;
    [pathlossmap.elevation_map pathlossmap.dtm_padding_ratio] = get_elevation_map_in_roi(pathlossmap, elevation_map_full, enable_debug_plotting);
    
   % if enable_debug_plotting
    if enable_debug_plotting
        coloraxis = [100,200];
        plot_pathloss_map(pathlossmap,coloraxis);
    end
    
    % Search for the index where to save this pathloss map (it is the
    % same order as how the eNodeBs are stored)
    filename_to_compare = strtok(strrep(pathlossmap.description.filename,'#','/'),'.');
    eNodeB_idx = [];
    sector_idx = [];
    for b_=1:length(eNodeBs)
        switch capesso_params.planning_tool
            case 'atoll'
                for s_=1:length(eNodeBs(b_).sectors)
                    site_name      = eNodeBs(b_).sectors(1).transmitter;
                    site_name_file = eNodeBs(b_).sectors(1).transmitter;
                    if strcmp(site_name,site_name_file)
                        eNodeB_idx = b_;
                        sector_idx = s_;
                        break
                    end
                end
            case 'capesso'
                site_name      = strtok(eNodeBs(b_).sectors(1).transmitter,'/');
                site_name_file = strtok(filename_to_compare,'/');
                if strcmp(site_name,site_name_file)
                    eNodeB_idx = b_;
                    sector_idx = 1;
                    break
                end
        end
        
    end
    
    % Save the pathloss file in the corresponding position
    if isempty(eNodeB_idx) || isempty(sector_idx)
        % Do not assign this pathloss to the output (no eNodeB/sector is using it)
        unused_pathloss_maps{length(unused_pathloss_maps)+1} = pathlossmap.description.filename;
    else
        % This site is not used
        if ~exist('pathlossmaps')
            pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
        else
            % This sector is not used
            if size(pathlossmaps,1)>=eNodeB_idx && size(pathlossmaps,2)>sector_idx
                if ~isempty(pathlossmaps(eNodeB_idx,sector_idx).description.filename)
                    unused_pathloss_maps{length(unused_pathloss_maps)+1} = pathlossmap.description.filename;
                else
                    pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
                end
            else
                pathlossmaps(eNodeB_idx,sector_idx) = pathlossmap;
             end
        end
    end
    
end

fclose('all');

% Assign the eNodeBs' sectors to the pathloss maps. Since these are isotrop maps, if
% we have maps for more than one sector (eg. W85A_32_t2#2FU1 and
% W85A_32_t2#2FU2), we will just take the #2FU (which means '/') and discard the rest.


function color_axis = plot_pathloss_map(pathlossmap_plot,varargin)
% In :
% pathlossmap_: containing data and description
% i_ for assignment to description data and subplot
%
% The extra argument allows you to specify the scale of the plot.
% Nevertheless, the used scale is returned by the function so you can
% use the same for all of your plots

figure_steps_col = 5;
figure_steps_row = 5;
pncols_ = pathlossmap_plot.description.ncols;
pnrows_ = pathlossmap_plot.description.nrows;
pxdim_ = pathlossmap_plot.description.xdim;
pydim_ = pathlossmap_plot.description.ydim;

% Attention : 3 and 4 are chosen static
% change, if more plots required

figure_capesso_los_file = figure('Colormap',[1 0 0;1 0.06275 0;1 0.1294 0;1 0.1922 0;1 0.2588 0;1 0.3216 0;1 0.3882 0;1 0.451 0;1 0.5176 0;1 0.5804 0;1 0.6471 0;1 0.7098 0;1 0.7725 0;1 0.8392 0;1 0.902 0;1 0.9686 0;0.9686 1 0;0.902 1 0;0.8392 1 0;0.7725 1 0;0.7098 1 0;0.6471 1 0;0.5804 1 0;0.5176 1 0;0.451 1 0;0.3882 1 0;0.3216 1 0;0.2588 1 0;0.1922 1 0;0.1294 1 0;0.06275 1 0;0 1 0;0 1 0.06275;0 1 0.1255;0 1 0.1882;0 1 0.251;0 1 0.3137;0 1 0.3765;0 1 0.4392;0 1 0.502;0 1 0.5608;0 1 0.6235;0 1 0.6863;0 1 0.749;0 1 0.8118;0 1 0.8745;0 1 0.9373;0 1 1;0 0.9373 1;0 0.8745 1;0 0.8118 1;0 0.749 1;0 0.6863 1;0 0.6235 1;0 0.5608 1;0 0.498 1;0 0.4392 1;0 0.3765 1;0 0.3137 1;0 0.251 1;0 0.1882 1;0 0.1255 1;0 0.06275 1;0 0 1]);
subplot(1,1, 1,'Parent',figure_capesso_los_file);

imagesc(pathlossmap_plot.pathloss_map(:,:));
set(gca,'YDir','normal');
set(gca,'XTick', 0 : pncols_/figure_steps_col: pncols_);
set(gca,'XTickLabel', 0:pncols_*pxdim_/figure_steps_col:pncols_*pxdim_);
set(gca,'YTick', 0 : pnrows_/figure_steps_row: pnrows_);
set(gca,'YTickLabel', 0:pnrows_*pydim_/figure_steps_row:pnrows_*pydim_);
xlabel('x pos [m]');
ylabel('y pos [m]');
if length(varargin)>=1
    color_axis = varargin{1};
    caxis(color_axis);
else
    color_axis = caxis;
end
colorbar;
title(['Capesso .LOS File: ',pathlossmap_plot.description.filename]);



function pathloss_files_description = get_pathloss_files_description(los_file_name, folder)
% Description of .LOS stored in .PAR Files

par_files = dir(fullfile(folder,'*.par'));
par_index = 0;

% Search for corresponding .PAR File in Folder
for j_=1:1:length(par_files)
    if strcmp(strtok(los_file_name,'.'),strtok(par_files(j_).name,'.'))
        par_index = j_;
        break;
    end
end

% Read out corresponding .PAR file
if(par_index)
    par_fid = fopen(fullfile(folder,par_files(par_index).name),'r');
    c = textscan(par_fid, '%s %f', 'delimiter','=');
    
    descr.filename = los_file_name;
    descr.NWxmap   = c{2}(1);
    descr.NWymap   = c{2}(2);
    descr.nrows    = c{2}(3);
    descr.ncols    = c{2}(4);
    descr.xdim     = c{2}(5);
    descr.ydim     = c{2}(6);
    descr.masked   = c{2}(7);
    
    descr.SWxmap = descr.NWxmap;                            % Lower-leftmost corner (x) -> SW
    descr.SWymap = descr.NWymap - descr.ydim*(descr.nrows-1/100); % Lower-leftmost corner (y) -> SW
    
    descr.NExmap = descr.NWxmap + descr.xdim*(descr.ncols-1/100);
    descr.NEymap = descr.NWymap;
    
    descr.SExmap = descr.NExmap;
    descr.SEymap = descr.SWymap;
    
    descr.roi_x = [descr.SWxmap descr.SExmap];
    descr.roi_y = [descr.SWymap descr.NWymap];
    
    fclose(par_fid);
    
    pathloss_files_description = descr;
else
    % TODO ... Return value
    % Implement print_log
    print_log(1,['Could not find .par file for ', los_file_name]);
end


function [aligned_elevation_map padding_ratio_p] = get_elevation_map_in_roi(pathlossmap, elevation_map_full,enable_debug_plotting)
% Returns the elevation map in the region of interest of the current pathlossmap.

pm_description  = pathlossmap.description;
pm_roi_min      = [pm_description.roi_x(1) pm_description.roi_y(1)];
pm_roi_max      = [pm_description.roi_x(2) pm_description.roi_y(2)];
dtm_description = elevation_map_full.description;
dtm_roi_min     = [elevation_map_full.description.roi_x(1) elevation_map_full.description.roi_y(1)];
dtm_roi_max     = [elevation_map_full.description.roi_x(2) elevation_map_full.description.roi_y(2)];

% Convert DTM coordinates (pixels) to PM coordinates (pixels)
pm_roi_min_pix_from_dtm = LTE_common_pos_to_pixel(pm_roi_min,dtm_roi_min,elevation_map_full.description.xdim);
pm_roi_max_pix_from_dtm = LTE_common_pos_to_pixel(pm_roi_max,dtm_roi_min,elevation_map_full.description.xdim);
dtm_roi_min_pix_from_pm = [1 1]; % Minimum pixel coordinates from the DTM map
dtm_roi_max_pix_from_pm = LTE_common_pos_to_pixel(dtm_roi_max,dtm_roi_min,elevation_map_full.description.xdim); % Maximum pixel coordinates from the DTM map

% The DTM should be bigger than the ROI, but it may not be the case
dtm_roi_min_pix_to_use = max(pm_roi_min_pix_from_dtm,dtm_roi_min_pix_from_pm);
dtm_roi_max_pix_to_use = min(pm_roi_max_pix_from_dtm,dtm_roi_max_pix_from_pm);

% Convert this square area to PM coordinates
dtm_in_pm_coord_min_pos = LTE_common_pixel_to_pos(dtm_roi_min_pix_to_use,dtm_roi_min,elevation_map_full.description.xdim);
dtm_in_pm_coord_max_pos = LTE_common_pixel_to_pos(dtm_roi_max_pix_to_use,dtm_roi_min,elevation_map_full.description.xdim);
dtm_in_pm_coord_min_pix = LTE_common_pos_to_pixel(dtm_in_pm_coord_min_pos,pm_roi_min,pm_description.xdim);
dtm_in_pm_coord_max_pix = LTE_common_pos_to_pixel(dtm_in_pm_coord_max_pos,pm_roi_min,pm_description.xdim);
dtm_cut_x_in_pm  = dtm_in_pm_coord_min_pix(1):dtm_in_pm_coord_max_pix(1);
dtm_cut_y_in_pm  = dtm_in_pm_coord_min_pix(2):dtm_in_pm_coord_max_pix(2);
dtm_cut_x_in_dtm = dtm_roi_min_pix_to_use(1):dtm_roi_max_pix_to_use(1);
dtm_cut_y_in_dtm = dtm_roi_min_pix_to_use(2):dtm_roi_max_pix_to_use(2);

% Put elevation map into the pm-aligned map
elevation_map_cut      = elevation_map_full.data(dtm_cut_y_in_dtm,dtm_cut_x_in_dtm);
elevation_map_cut_mean = mean(elevation_map_cut(:));
aligned_elevation_map  = zeros(pathlossmap.description.nrows, pathlossmap.description.ncols)+elevation_map_cut_mean; % Padding is done with the average value
aligned_elevation_map(dtm_cut_y_in_pm,dtm_cut_x_in_pm) = elevation_map_cut;

% Calculation of what was put
elements_put    = numel(elevation_map_cut);
total_elements  = numel(aligned_elevation_map);
padded_elements = total_elements-elements_put;
padding_ratio_p = padded_elements/total_elements * 100;

% if padded_elements
%     warning('DTM did not cover the whole pathloss map, padding was added. %3.2f%% of the area corresponding to the pathloss map was padded at %3.1fm\n',padding_ratio_p,elevation_map_cut_mean);
% end

% figure;imagesc(pm_description.roi_x,pm_description.roi_y,cut_elevation_map);

if enable_debug_plotting
    % Plot whole elevation map
    figure;
    imagesc(elevation_map_full.description.roi_x,elevation_map_full.description.roi_y,elevation_map_full.data);
    color_axis = caxis;
    grid on;
    hold on;
    set(gca,'YDir','normal');
    xlabel('x pos [m]');
    ylabel('y pos [m]');
    % Mark the part of the elevation map that corresponds to the pathloss file
    roi_width  = pm_description.roi_x(2)-pm_description.roi_x(1);
    roi_height = pm_description.roi_y(2)-pm_description.roi_y(1);
    rectangle('Position',[pm_roi_min(1), pm_roi_min(2), pm_roi_max(1)-pm_roi_min(1), pm_roi_max(2)-pm_roi_min(2)]);
    hold off;
    colorbar;
    title('Elevation map (m)');
    % Plot only area corresponding to this pathloss map
    figure;
    imagesc([pm_roi_min(1) pm_roi_max(1)],[pm_roi_min(2) pm_roi_max(2)],aligned_elevation_map);
    caxis(color_axis);
    grid on;
    set(gca,'YDir','normal');
    colorbar;
    title(sprintf('Elevation map (m). %s area',strrep(pm_description.filename,'_','\_')));
    xlabel('x pos [m]');
    ylabel('y pos [m]');
end


