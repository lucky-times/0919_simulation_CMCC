function [eNodeBs eNodeBs_sectors networkPathlossMap networkShadowFadingMap] = LTE_init_network_generation
%% Network generation. Either from a file (cache) or calling the necessary function.
% (c) Josep Colom Ikuno, INTHFT, 2009
% www.nt.tuwien.ac.at

global LTE_config;

% If cache is on and the cache file exists, load network from disk
if LTE_config.cache_network && exist(LTE_config.network_cache,'file')
    print_log(1,sprintf('Loading network from %s\n',LTE_config.network_cache));
    load(LTE_config.network_cache);
    % Change map resolution to match the loaded maps'
    LTE_config.map_resolution = networkPathlossMap.data_res;
% If not, then create it
else 
    % Generate network (eNodeBs and macroscopic pathloss)
    print_log(1,'Generating network\n');
    switch LTE_config.network_source
        case 'generated'
            [eNodeBs networkPathlossMap] = LTE_init_generate_network;
            % Other case here -> other sources, eg. Network planning tool
        case 'capesso'
            [eNodeBs networkPathlossMap] = LTE_init_generate_capesso_network;
            display('Using data from planning tool Capesso');
        otherwise
            error([LTE_config.network_source ' network source not supported\n']);
    end
    
    % Add the power separation. X% to signaling/pilots (always on) and the rest for data
    for b_=1:length(eNodeBs)
        for s_=1:length(eNodeBs(b_).sectors)
            data_power      = eNodeBs(b_).sectors(s_).max_power * (1-LTE_config.signaling_ratio);
            signaling_power = eNodeBs(b_).sectors(s_).max_power * LTE_config.signaling_ratio;
            eNodeBs(b_).sectors(s_).max_power       = data_power;
            eNodeBs(b_).sectors(s_).signaling_power = signaling_power;
            LTE_config.scheduler_params.max_power   = data_power; % max data transmit power in Watts
        end
    end
    
    % A posteriori calculation of neighboring eNodeBs
    for b_ = 1:length(eNodeBs)
        % Use the 6 closest sites as interferers
        [eNodeBs(b_).neighbors eNodeBs(b_).neighbors_eNodeB] = LTE_init_get_eNodeB_neighbors(eNodeBs(b_), eNodeBs);
    end
    
    % Generate shadow fading
    if strcmp(LTE_config.network_source,'generated')
        print_log(1,'Generating shadow fading\n');
        switch LTE_config.shadow_fading_type
            case 'claussen'
                [LTE_config.roi_x LTE_config.roi_y] = networkPathlossMap.valid_range;
                networkShadowFadingMap = LTE_init_generate_claussen_shadow_fading_map(eNodeBs);
            case 'none'
                [LTE_config.roi_x LTE_config.roi_y] = networkPathlossMap.valid_range;
                networkShadowFadingMap = LTE_init_generate_NoShadowing_shadow_fading_map(eNodeBs);
            otherwise
                error([LTE_config.shadow_fading_type ' shadow fading type not supported']);
        end
    else
      %  error('only "generated" supported for now');
    end
    
    %% Fill in sector assignment (takes into account the shadow fading)
    if exist('networkShadowFadingMap','var')
        networkPathlossMap.sector_assignment = LTE_common_calculate_sector_assignment(networkPathlossMap);
        networkPathlossMap.sector_assignment = LTE_common_calculate_sector_assignment(networkPathlossMap,networkShadowFadingMap);
    else
        networkPathlossMap.sector_assignment = LTE_common_calculate_sector_assignment(networkPathlossMap);
    end
    
    % Save network
    if LTE_config.cache_network
        print_log(1,'Saving network to file\n');
        if exist('networkShadowFadingMap','var')
            save(LTE_config.network_cache,'eNodeBs','networkPathlossMap','networkShadowFadingMap');
        else
            save(LTE_config.network_cache,'eNodeBs','networkPathlossMap');
        end
    end
end

% Reapply MCL to the pathloss maps (not saved in the cache)
if LTE_config.minimum_coupling_loss
    networkPathlossMap.pathloss = max(networkPathlossMap.pathloss,LTE_config.minimum_coupling_loss);
end

% Calculate eNodeBs' distances
eNodeBs_positions = reshape([eNodeBs.pos],length(eNodeBs),[]);
distances = zeros(length(eNodeBs));
for b_=1:length(eNodeBs)
    distance = sqrt((eNodeBs_positions(b_,1)-eNodeBs_positions(:,1)).^2 + (eNodeBs_positions(b_,2)-eNodeBs_positions(:,2)).^2);
    distances(:,b_) = distance;
end
sorted_distances = sort(distances);
if length(eNodeBs) > 1
    max_inter_eNodeB_distance  = max(sorted_distances(2,:));
    min_inter_eNodeB_distance  = min(sorted_distances(2,:));
    mean_inter_eNodeB_distance = mean(sorted_distances(2,:));
else 
    max_inter_eNodeB_distance  = 0;
    min_inter_eNodeB_distance  = 0;
    mean_inter_eNodeB_distance = 0;
end

% Add number of antennas information, as well as sector info
eNodeBs_sectors = eNodeBs(1).sectors(1); % Initialization
eNodeBs_idx = 1;
for b_=1:length(eNodeBs)
    for s_=1:length(eNodeBs(b_).sectors)
        eNodeBs(b_).sectors(s_).nTX = LTE_config.nTX;
        eNodeBs(b_).sectors(s_).eNodeB_id = eNodeBs_idx;
        eNodeBs_sectors(eNodeBs_idx) = eNodeBs(b_).sectors(s_);
        eNodeBs_idx = eNodeBs_idx + 1;
    end
end

% Configure the case for zero-delay
if LTE_config.feedback_channel_delay==0
    for b_=1:length(eNodeBs)
        for s_=1:length(eNodeBs(b_).sectors)
            eNodeBs(b_).sectors(s_).zero_delay_feedback = true;
        end
    end
end

% Configure unquantized feedback
if LTE_config.unquantized_CQI_feedback
    for b_=1:length(eNodeBs)
        for s_=1:length(eNodeBs(1).sectors)
            eNodeBs(b_).sectors(s_).unquantized_CQI_feedback = true;
        end
    end
end

if ~exist('networkShadowFadingMap','var')
    % To avoid error of missing return argument
    networkShadowFadingMap = 0;
end
