function M_Capesso_return = LTE_get_antenna_gain_map(capesso_params,eNodeBs,M_Capesso)
% Generates the antenna gain map. Projects the interpolated 3D antenna pattern onto the elevation map.
% The Antenna directional radio pattern is read from a MSI file (eg. from Kathrein antenna files).
%
% (c) Martin Taranetz, Josep Colom Ikuno INTHFT, 2010

% Generate Matrix with antenna gain pattern for every eNodeB-Sector (1..3,:,:) and add it to M_Capesso

for b_ = 1:length(M_Capesso)
    % Geometrical data from pathloss map
    ulxmap_ = M_Capesso(b_).description.NWxmap;
    ulymap_ = M_Capesso(b_).description.NWymap;
    nrows_  = M_Capesso(b_).description.nrows;
    ncols_  = M_Capesso(b_).description.ncols;
    xdim_   = M_Capesso(b_).description.xdim;
    ydim_   = M_Capesso(b_).description.ydim;
    
    M_Capesso(b_).antenna_gain_map = zeros(nrows_, ncols_, length(eNodeBs(b_).sectors));
    
    eNodeB_pos      = eNodeBs(b_).pos;
    eNodeB_altitude = eNodeBs(b_).altitude;
    
    % Build two matrices for x positions and y positions
    % Allthough redundant this enables an effective implementation of
    % calculating the angle maps.
    position_grid_meters = zeros(nrows_,ncols_,2);
    position_grid_meters(:,:,1) = repmat((ulxmap_:xdim_:ulxmap_+xdim_*ncols_-1), nrows_, 1);
    
    % As terrain- and pathlossmaps are stored such that (1,1) refers to SW
    % point of the map, we have to generate the horizontal map accordingly
    position_grid_meters(:,:,2) = repmat((ulymap_-ydim_*(nrows_-1):ydim_:ulymap_)', 1, ncols_);
    
    % HORIZONTAL angle grid
    horizontal_angle_grid   = rad2deg(atan2((position_grid_meters(:,:,1) - eNodeB_pos(1)),(position_grid_meters(:,:,2)-eNodeB_pos(2))));
    horizontal_angle_grid_w = wrapTo360(horizontal_angle_grid);
    % horizontal distance to eNodeB
    distance_map = sqrt((position_grid_meters(:,:,1) - eNodeB_pos(1)).^2 + (position_grid_meters(:,:,2)-eNodeB_pos(2)).^2);
    % VERTICAL angle grid
    tx_heights = unique([eNodeBs(b_).sectors.tx_height]);
    vertical_angle_grid_el = zeros([size(M_Capesso(b_).elevation_map) length(tx_heights)]);
    if capesso_params.enable_dtm
        for tx_heights_idx=1:length(tx_heights)
            vertical_angle_grid_el(:,:,tx_heights_idx) = rad2deg(atan2((eNodeB_altitude + tx_heights(tx_heights_idx) - (M_Capesso(b_).elevation_map + capesso_params.rx_height)), distance_map));
        end
    else
        % Do not use DTM data. Flat scenario assumed - still considering
        % transmitter and receiver height
        for tx_heights_idx=1:length(tx_heights)
            vertical_angle_grid_el(:,:,tx_heights_idx) = rad2deg(atan2((tx_heights(tx_heights_idx) - capesso_params.rx_height), distance_map));
        end
    end

    if capesso_params.enable_debug_plotting
        figure;
    end
    % Generate gain pattern for every sector of eNodeB
    for s_ = 1:length(eNodeBs(b_).sectors)
        azimuth = eNodeBs(b_).sectors(s_).azimuth;
        horizontal_angle_grid_s = wrapTo360(horizontal_angle_grid_w - azimuth);
        % Use mechanical tilt values from Capesso data or default value -
        antenna_electrical_tilt   = eNodeBs(b_).sectors(s_).electrical_downtilt;
        if capesso_params.use_default_tilt_value == false
            antenna_mechanical_tilt   = eNodeBs(b_).sectors(s_).mechanical_downtilt;
        else
            antenna_mechanical_tilt   = capesso_params.default_mechanical_tilt;
        end    
        vertical_angle_grid_el_idx = find(tx_heights==eNodeBs(b_).sectors(s_).tx_height); % unique values. No need to take duplicates into acocunt
        gain_pattern_ = eNodeBs(b_).sectors(s_).antenna.gain(horizontal_angle_grid_s, vertical_angle_grid_el(:,:,vertical_angle_grid_el_idx), antenna_electrical_tilt, antenna_mechanical_tilt);
          
        M_Capesso(b_).antenna_gain_map(:,:,s_) = gain_pattern_;
        % Plot antenna gain maps
        if capesso_params.enable_debug_plotting
            subplot(2,2,s_);
            imagesc(gain_pattern_);
            colorbar;
            title(strrep(sprintf('%s (Sec%d), \\theta=%3.0f°',strrep(M_Capesso(b_).description.filename,'#2','/'),s_,azimuth),'_','\_'));
            caxis([-60,20]);
        end
    end
end

M_Capesso_return = M_Capesso;

