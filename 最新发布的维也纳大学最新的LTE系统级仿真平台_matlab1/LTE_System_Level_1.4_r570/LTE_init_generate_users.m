function [ UEs extra_info ] = LTE_init_generate_users(eNodeBs,networkMacroscopicPathlossMap)
% Get a cell assignment matrix. Here I will directly take it from the
% networkMacroscopicPathlossMap, but the following code could generate an equivalent
% one for any type of pathloss model.
%
% global LTE_config;
% [ x_range y_range ] = networkMacroscopicPathlossMap.valid_range;
% step = LTE_config.map_resolution;
% array_x =min(x_range):step:max(x_range);
% array_y =min(y_range):step:max(y_range);
% loss = zeros(length(array_y),length(array_x));
% for x_ = 1:length(array_x)
%     for y_ = 1:length(array_y)
%         pos = [array_x(x_) array_y(y_)];
%         assignment_map(y_,x_) = networkMacroscopicPathlossMap.cell_assignment(pos(1),pos(2));
%     end
% end
%
% (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010

global LTE_config;

use_UE_cache = LTE_config.UE_cache;
cache_file_exists = exist(LTE_config.UE_cache_file,'file');


%% Data needed also for plotting: generate every time

if (LTE_config.UE.use_traffic_map==false)

    sector_surfaces = networkMacroscopicPathlossMap.sector_sizes;
    
    % Check whether some sector surfaces are zero. This could be caused by
    % choosing a very small inter-eNodeB distance in relation to the map
    % resolution.
    if sum(sector_surfaces(:)==0)>0
        warning('Some sector sizes are zero. No UEs will be generated there. Maybe a too big map resolution value?');
    end
    
    % Calculate how many users per sector
    % norm_sector_surface = sector_surfaces / sector_surfaces(LTE_config.target_sector(1),LTE_config.target_sector(2));
    % Before the UEs per sector were normalized to the size of the target sector. This is now changed to have a fixed amount per sector
    norm_sector_surface = ones(size(sector_surfaces));
    users_sector = round(norm_sector_surface*LTE_config.UE_per_eNodeB).*(sector_surfaces>0);
    % Alternatively use UE_per_eNodeB as threshold probability that single
    % user is generated in sector: UE_per_eNodeB in [0 1]
    % user_sector = rand(size(norm_sector_surface))<UE_per_eNodeB;
    
    % All the possible positions for a given sector
    sector_positions = cell(size(networkMacroscopicPathlossMap.sector_sizes));
    % Where our users will be positions (for each eNodeB and sector)
    user_positions_pixels = cell(size(networkMacroscopicPathlossMap.sector_sizes));
    % Assign random positions to each UE
    sector = 1;
    bts    = 2;
    for b_ = 1:length(eNodeBs)
        for s_ = 1:length(eNodeBs(1).sectors)
            if users_sector(b_,s_)~=0
                sector_positions_matrix = networkMacroscopicPathlossMap.sector_assignment(:,:,sector)==s_ & networkMacroscopicPathlossMap.sector_assignment(:,:,bts)==b_;
                [row,col] = find(sector_positions_matrix);
                sector_positions{b_,s_} = [col,row];
                user_positions_pixels{b_,s_} = sector_positions{b_,s_}(ceil(size(sector_positions{b_,s_},1)*rand(1,users_sector(b_,s_))),:);
            else
                sector_positions{b_,s_} = [];
                user_positions_pixels{b_,s_} = [];
            end
        end
    end
else
    % Use traffic maps for user generation
    
    % OUR TESTCASE :  
    % user density per pixel
    % user_density_ = 0.02;
    % traffic_map_in_roi = repmat(user_density_, size(networkMacroscopicPathlossMap.sector_assignment(:,:,1)));
    
    % Using User Density Traffic Maps from MKA
    % Last parameter enables or disables plotting
    user_density_traffic_maps = LTE_init_user_density_traffic_maps(LTE_config.udtm_folder, LTE_config.udtm_filename, LTE_config.udtm_environment, false);
    % Set undefined values to zero
    user_density_traffic_maps.data(user_density_traffic_maps.data<0) = 0;
    % Now rescale the user density traffic map and calculate a user
    % density/pixel map
    
    % Resize by rescale factor also used for capesso pathlossmaps
    user_density_traffic_maps.data = imresize(user_density_traffic_maps.data , LTE_config.rescale_factor);
    % Cut out the ROI
    pm.roi_x = networkMacroscopicPathlossMap.roi_x;
    pm.roi_y = networkMacroscopicPathlossMap.roi_y;
    udtm.roi_x = user_density_traffic_maps.description.roi_x;
    udtm.roi_y = user_density_traffic_maps.description.roi_y;
    % Calculate SW position of Pathlossmaps in pixel within Traffic maps
    pm_pos = LTE_common_pos_to_pixel([pm.roi_x(:,1) pm.roi_y(:,1)],[udtm.roi_x(:,1) udtm.roi_y(:,1)], networkMacroscopicPathlossMap.data_res);
    n_rows = size(networkMacroscopicPathlossMap.pathloss,1);
    n_cols = size(networkMacroscopicPathlossMap.pathloss,2);
    traffic_map_in_roi_pkm = user_density_traffic_maps.data(pm_pos(2):pm_pos(2)+n_rows-1, pm_pos(1):pm_pos(1)+n_cols-1);
    
    % Recalculate the traffic map in the ROI to a users/pixel map
    % Consider both dimensions to be equal:
    x_dim_n = user_density_traffic_maps.description.xdim / LTE_config.rescale_factor;
    square_meters_per_pixel = x_dim_n^2;
    traffic_map_in_roi_pp = traffic_map_in_roi_pkm/(10^6) * square_meters_per_pixel;
    
    % Place users in ROI via coin toss - use users/pixel density as
    % probability for positive outcome of random experiment

    % using traffic_map_in_roi as threshold and uniform distribution gives
    % logical matrix of pixel user positions in ROI
    user_matrix_in_roi = (rand(size(traffic_map_in_roi_pp)) <= traffic_map_in_roi_pp);
    
    % Find eNodeB and sector for each generated user and put position into
    % matrix element according to sector_assignment
    
     % Where our users will be positions (for each eNodeB and sector)
    user_positions_pixels = cell(size(networkMacroscopicPathlossMap.sector_sizes));
    % All the possible positions for a given sector
    sector_positions = cell(size(networkMacroscopicPathlossMap.sector_sizes));
    
    sector = 1;
    bts    = 2;
    for b_ = 1:length(eNodeBs)
        for s_ = 1:length(eNodeBs(1).sectors)
            % Used to mask out users of appropriate sector and eNodeB
            sector_positions_matrix = networkMacroscopicPathlossMap.sector_assignment(:,:,sector)==s_ & networkMacroscopicPathlossMap.sector_assignment(:,:,bts)==b_;
            [row, col] = find(sector_positions_matrix.*user_matrix_in_roi);
            [row_sp,col_sp] = find(sector_positions_matrix);
            sector_positions{b_,s_} = [col_sp, row_sp];
            if ~isempty(row)
                user_positions_pixels{b_,s_} = [col, row];
            else
                user_positions_pixels{b_,s_} = [];
            end
        end
    end
end

%% Creating or loading UE position, depending on the configuration
% Create UEs according to the previously generated positions
UEs = network_elements.UE;
if (~use_UE_cache) || (use_UE_cache&&~cache_file_exists)    
    userID = 1;
    for b_ = 1:length(eNodeBs)
        for s_ = 1:length(eNodeBs(1).sectors)
            
            % Generate only necessary users
            if ~LTE_config.UEs_only_in_target_sector
                UEs_to_create = 1:size(user_positions_pixels{b_,s_},1);
            else
                if b_==LTE_config.target_sector(1) && s_==LTE_config.target_sector(2)
                    UEs_to_create = 1:size(user_positions_pixels{b_,s_},1);
                else
                    UEs_to_create = [];
                end
            end
            
            UE_positions = zeros(userID,2);
            for u_ = UEs_to_create
                % General UE settings that can be saved and re-used
                UEs(userID)     = network_elements.UE;
                UEs(userID).id  = userID;
                UEs(userID).pos = LTE_common_pixel_to_pos( user_positions_pixels{b_,s_}(u_,:), networkMacroscopicPathlossMap.coordinate_origin, networkMacroscopicPathlossMap.data_res);
                % Generate a walking model for the user
                if isfield(LTE_config.UE,'walk')
                    if strcmp(LTE_config.UE.walk,'SLvsLL')
                        UEs(userID).walking_model = walking_models.SLvsLLWalkingModel; % Since no angle is specified, a random one is chosen
                    end
                else
                    UEs(userID).walking_model = walking_models.straightWalkingModel(LTE_config.UE_speed*LTE_config.TTI_length); % Since no angle is specified, a random one is chosen
                end
                UE_positions(userID,:) = user_positions_pixels{b_,s_}(u_,:);
                if isfield(LTE_config.UE,'walk')
                    if strcmp(LTE_config.UE.walk,'SLvsLL')
                        UE_positions(userID,:) = [5,8.6];
                        temp_pos = LTE_common_pos_to_pixel(UE_positions(userID,:), networkMacroscopicPathlossMap.coordinate_origin, networkMacroscopicPathlossMap.data_res);
                        UEs(userID).pos = LTE_common_pixel_to_pos( temp_pos, networkMacroscopicPathlossMap.coordinate_origin, networkMacroscopicPathlossMap.data_res);
                    end
                end
                
                userID = userID + 1;
            end
        end
    end
    print_log(1,sprintf('Saving UE positions to %s\n',LTE_config.UE_cache_file));
    save(LTE_config.UE_cache_file,'UEs','UE_positions');
else
    % Load UEs
    print_log(1,sprintf('Loading UE positions from %s\n',LTE_config.UE_cache_file));
    UE_data_from_cache = load(LTE_config.UE_cache_file);
    UEs = UE_data_from_cache.UEs;
    
    % Get a struct name with everything except the UEs field
    names = fieldnames(UE_data_from_cache);
    if length(names) > 1 % If there is nothing more besides the 'UEs' data
        extra_info = struct;
        for field_idx=1:length(names)
            current_name = names{field_idx};
            if ~strcmp(current_name,'UEs')
                extra_info = setfield(extra_info,current_name,getfield(UE_data_from_cache,current_name));
            end
        end
    end
end

% Safeguard against having no UEs
if length(UEs)==1 && isempty(UEs(1).id)
    no_UEs = true;
else
    no_UEs = false;
end

if no_UEs
    output_args = [];
else
    % Assign the UEs to their nearest (in pathloss) eNodeB and assign some extra parameters
    UE_positions_m = zeros(length(UEs),2);
    for u_ = 1:length(UEs)
        % To be sure: assign each user according to the cell assignment map
        [ eNodeB_id sector_num ] = networkMacroscopicPathlossMap.cell_assignment(UEs(u_).pos);
        UEs(u_).attached_eNodeB = eNodeBs(eNodeB_id);
        UEs(u_).attached_sector = sector_num;
        
        % Set the channel model for the user
        % This is now done outside of the function
        
        % Attach UE to eNodeB
        eNodeBs(eNodeB_id).attachUser(UEs(u_),sector_num);
        UE_positions_m(u_,:) = UEs(u_).pos;
        
        % Append traffic model to users
        UEs(u_).traffic_model = LTE_trafficmodel(LTE_config.traffic_models,UEs(u_),max(LTE_config.feedback_channel_delay,0));
        UEs(u_).trace_SINR    = false; %LTE_config.traces_config.trace_SINR;
    end
    
    %% Plot where the users are (pixel positions)
    if LTE_config.show_network>0
        figure(LTE_config.plots.initial_UE_positions);
        hold on;
        roi_min = [min(networkMacroscopicPathlossMap.roi_x) min(networkMacroscopicPathlossMap.roi_y)];
        total_elements = length(eNodeBs)*length(eNodeBs(1).sectors);
        h = (1:total_elements)/total_elements;
        % the best way I could find to generate a randomly selected saturated colormap for the different sectors
        colormaps = hsv2rgb([h' ones(length(h),2)]);
        colormap_permutation = randperm(size(colormaps,1));
        colormaps = colormaps(colormap_permutation,:);
        i_=1;
        for b_ = 1:length(eNodeBs)
            for s_ = 1:length(eNodeBs(1).sectors)
                if ~isempty(sector_positions{b_,s_}) % An error would happen if no UEs are present
                    current_sector_positions = LTE_common_pixel_to_pos(sector_positions{b_,s_},roi_min,networkMacroscopicPathlossMap.data_res);
                    scatter(current_sector_positions(:,1),current_sector_positions(:,2),'+','MarkerEdgeColor',colormaps(i_,:),'MarkerFaceColor',colormaps(i_,:));
                end
                i_ = i_+1;
            end
        end
        scatter(UE_positions_m(:,1),UE_positions_m(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','w');
        
        title(sprintf('UE initial positions: %d eNodeBs, %d sectors/eNodeB ',length(eNodeBs),length(eNodeBs(1).sectors)));
        xlim(networkMacroscopicPathlossMap.roi_x);
        ylim(networkMacroscopicPathlossMap.roi_y);
        xlabel('x pos [m]');
        ylabel('y pos [m]');
    end
    
    % Choose as many points per cell as users
    if ~exist('extra_info','var')
        extra_info = [];
    end
end
