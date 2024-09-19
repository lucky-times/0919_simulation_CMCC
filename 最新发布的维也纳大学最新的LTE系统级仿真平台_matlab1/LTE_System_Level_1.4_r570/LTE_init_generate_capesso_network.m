function [eNodeBs networkMacroscopicPathlossMap] = LTE_init_generate_capesso_network
% (c) Martin Taranetz , INTHFT, 2010

% output:   eNodeBs            ... contains info reagarding the BTSs and its
%                                  sectors
%           pathloss_data      ... [heightxwidthx3xnBTS] double
%                                  Pathloss data for each sector (including
%                                  antenna gain).
%                                  [y,x,sector_num,brts_num]

%% Configuration parameters

% Necessary data from the LTE_config
global LTE_config;

capesso_params = LTE_config.capesso_params;
capesso_params.frequency                  = LTE_config.frequency;
capesso_params.maps_resolution            = LTE_config.map_resolution; 
capesso_params.rescale_factor             = LTE_config.rescale_factor;
capesso_params.enable_plotting            = LTE_config.show_network;

%% Create the eNodeBs
eNodeBs   = LTE_init_read_cell_data(capesso_params);

%% Create the Pathlossmaps

print_log(1,'Creating cell pathloss map from Capesso data\n');
M_Capesso = LTE_get_capesso_pathlossmaps(capesso_params, eNodeBs); % Since this data is from an isotrop antenna, just one map per site is loaded.
print_log(1,'Calculating sector antenna gains\n');
M_Capesso = LTE_get_antenna_gain_map(capesso_params, eNodeBs, M_Capesso); % Add antenna gain

networkMacroscopicPathlossMap                        = channel_gain_wrappers.macroscopicPathlossMap;
networkMacroscopicPathlossMap.data_res               = capesso_params.maps_resolution / capesso_params.rescale_factor ;

%% Calculate ROI and find BTS which is nearest to the middle of the ROI
% BTS positions
tx_pos = zeros(length(eNodeBs),2);
for b_ = 1:length(eNodeBs)
    tx_pos(b_,:) = eNodeBs(b_).pos;
end

%% Manually set ROI for drive tests
if ~LTE_config.manually_set_ROI
    print_log(1,'ROI bounded by the eNodeBs plus increase factor');
    % Calculate ROI border points in ABSOLUTE coordinates
    roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
    roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
    % roi_reduction_factor times smaller and draw it. ABSOLUTE COORDINATES
    roi_x = roi_x + capesso_params.eNodeB_ROI_increase_factor*abs(roi_x(2)-roi_x(1))*[-1,1];
    roi_y = roi_y + capesso_params.eNodeB_ROI_increase_factor*abs(roi_y(2)-roi_y(1))*[-1,1];
else
    print_log(1,'ROI set manually');
    roi_x = LTE_config.roi_x;
    roi_y = LTE_config.roi_y;
end

networkMacroscopicPathlossMap.roi_x = roi_x;
networkMacroscopicPathlossMap.roi_y = roi_y;

%%

% Find the BTS that is nearest to the middle of the ROI
if strcmp(LTE_config.target_sector,'center')
    middle_ = [mean(roi_x),mean(roi_y)];
    target_bts = 1;
    for i_ = 1:length(eNodeBs)
        if norm(tx_pos(i_,:)-middle_) < norm(tx_pos(target_bts,:)-middle_)
            target_bts = i_;
        end
    end
    LTE_config.target_sector    = target_bts;
    LTE_config.target_sector(2) = 1;
else
    target_bts = LTE_config.target_sector(1);
end
eNodeBs(target_bts).target_cell = true;

%% Calculate sector pathloss in ROI

% Size of ROI in pixels
roi_size_in_pixels = LTE_common_pos_to_pixel([roi_x(:,2) roi_y(:,2)],[roi_x(:,1) roi_y(:,1)], capesso_params.maps_resolution);
roi_width_pixels = roi_size_in_pixels(:,1);
roi_height_pixels = roi_size_in_pixels(:,2);

% Memory preallocation
sector_pathloss_data = zeros(roi_height_pixels*capesso_params.rescale_factor , roi_width_pixels*capesso_params.rescale_factor, capesso_params.number_of_sectors, length(eNodeBs));

% Calculate final pathloss in the region of interest for every sector
% applying the antenna gain
print_log(1,'Creating sector pathloss map (applying sector antenna gains)\n');

for b_=1:length(eNodeBs)
   % Calculate pixel position of ROI in pathlossmap
   % Position in pixel is the same for pathlossmap and antenna gain
   % patterns of sectors
   roi_position_in_pixels = LTE_common_pos_to_pixel([roi_x(:,1) roi_y(:,1)], [M_Capesso(b_).description.SWxmap M_Capesso(b_).description.SWymap]  , capesso_params.maps_resolution);
   % Get pathlossmap in region of interest
   pathlossmap_in_roi = M_Capesso(b_).pathloss_map(roi_position_in_pixels(2):roi_position_in_pixels(2)+roi_height_pixels-1, roi_position_in_pixels(1):roi_position_in_pixels(1)+roi_width_pixels-1);
   pathlossmap_in_roi = imresize(pathlossmap_in_roi, capesso_params.rescale_factor);
   
   % Store pathlosses
   if b_==1
       pathlossmap_in_roi_all = zeros([size(pathlossmap_in_roi), length(eNodeBs)]);
   end
   pathlossmap_in_roi_all(:,:,b_) = pathlossmap_in_roi;
   
   % Debug plot of pathlossmap in ROI
   if capesso_params.enable_debug_plotting
     figure;
     imagesc(pathlossmap_in_roi);
     title([M_Capesso(b_).description.filename num2str(b_)]);
     set(gca,'YDir','normal');
     colorbar;
   end

    for s_=1:length(eNodeBs(b_).sectors)
       % Antenna gain in region of interest
       sector_antenna_gain_in_roi = M_Capesso(b_).antenna_gain_map(roi_position_in_pixels(2):roi_position_in_pixels(2)+roi_height_pixels-1, roi_position_in_pixels(1):roi_position_in_pixels(1)+roi_width_pixels-1, s_);
       sector_antenna_gain_in_roi = imresize(sector_antenna_gain_in_roi, capesso_params.rescale_factor);
       % Calculate final pathloss without minimum coupling loss
       sector_pathloss_data(:,:,s_,b_) = pathlossmap_in_roi - sector_antenna_gain_in_roi;
       
%      Display min/max values for debugging
%        print_log(1,min(min(pathlossmap_in_roi)));
%        print_log(1,max(max(pathlossmap_in_roi)));
%        print_log(1,min(min(sector_antenna_gain_in_roi)));
%        print_log(1,max(max(sector_antenna_gain_in_roi)));
%        print_log(1,min(min(sector_pathloss_data(:,:,s_,b_))));
%        print_log(1,max(max(sector_pathloss_data(:,:,s_,b_))));
       
   end
end

% Plot isotropic SIR (debug only)
if capesso_params.enable_debug_plotting
    plot_capesso_isotropic_SINR(pathlossmap_in_roi_all);
    plot_capesso_isotropic_SINR(sector_pathloss_data);
end

% Fill in pathlossdata in the pathloss map
networkMacroscopicPathlossMap.pathloss = sector_pathloss_data;

%% Plotting
% Creates Figures for each eNodeB containing
%   o  Capesso pathloss map
%   o  Elevation map
%   o  eNodeB positions
%   o  Antenna gain map for each sector
if capesso_params.enable_plotting
    LTE_plot_capesso_files(eNodeBs, M_Capesso, networkMacroscopicPathlossMap, capesso_params);
end

% Debug function to plot the approximate shape of the cells
function plot_capesso_isotropic_SINR(isotropic_pathloss)
% Convert to 3D
if ndims(isotropic_pathloss) > 3
    isotropic_pathloss = reshape(isotropic_pathloss,size(isotropic_pathloss,1),size(isotropic_pathloss,2),[]);
end

% Calculate number of files
n_files = size(isotropic_pathloss,3);

pathloss_lin = 10.^(-isotropic_pathloss/10);
for b_=1:n_files
    others = [1:b_-1 b_+1:n_files];
    SIR_lin(:,:,b_) = pathloss_lin(:,:,b_) ./ sum(pathloss_lin(:,:,others),3);
end
SIR_dB = 10*log10(SIR_lin);

for b_=1:n_files
    figure;
    climits = [min(SIR_dB(:)) max(SIR_dB(:))];
    imagesc(SIR_dB(:,:,b_),climits);
    set(gca,'YDir','normal');
    colorbar;
    title(sprintf('Site SIR [dB], clipped to [%3.1f %3.1f] dB',climits(1),climits(2)));
end

[max_SIR I] = max(SIR_dB,[],3);
figure;imagesc(max_SIR); colorbar; set(gca,'YDir','normal'); title('Isotropic maximum SIR [dB]');
figure;imagesc(I); colorbar; set(gca,'YDir','normal'); title('isotropic cell allocation [dB]');