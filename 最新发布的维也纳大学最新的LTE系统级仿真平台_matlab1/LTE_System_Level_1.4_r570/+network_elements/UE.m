classdef UE < handle
    % Class that represents an LTE UE (user)
    % (c) Josep Colom Ikuno, INTHFT, 2008
    
    properties
        % Unique UE id
        id
        % pos in meters (x,y)
        pos
        % eNodeB site to where this user ir attached
        attached_eNodeB
        % eNodeB index to which the UE is attached
        attached_sector
        % Walking model for this user
        walking_model
        % Downlink channel model for this user
        downlink_channel
        RB_grid % this links to obj.downlink_channel.RB_grid Stored again for efficiency reasons
        
        % Uplink channel model for this user
        uplink_channel
        % noise figure for this specific UE
        receiver_noise_figure
        % penetration loss in dB for this user
        penetration_loss
        % Number of receive antennas
        nRX
        % Antenna gain of the UE
        antenna_gain
        
        % trace that stores info about what happened
        trace
        % network clock. Tells the UE in what TTI he is
        clock
        % performs the mapping between SINR and CQI
        CQI_mapper
        
        % Output of the link quality (measurement) model. Contains the
        % following data:
        %  - nCodewords: how many codewords were used
        %  - CQI_feedback: the measured CQI
        link_quality_model_output % SINR_dB
                                  % SINR_linear

        % Data to be fed back to the eNodeB. It is used to pass the feedback data to the send_feedback() function
        feedback                  % tx_mode -> to differentiate when there is RI feedback
                                  % nCodewords
                                  % CQI
                                  % RI
                                  % TB_ACK
                                  % TB_size
                                  % TB_was_scheduled
        
        % Whether the CQI feedback should be unquantized. Having this set
        % to true is equivalent to directly sending the SINR
        unquantized_CQI_feedback = false;
        
        % Will decide whether a give TB made it or not
        BLER_curves
        
        % Gives the means to average the several Transport Block (TB) SINRs
        SINR_averager
        
        % Signaling from the eNodeB to this UE. This is a direct channel
        % between the eNodeB and this UE, where it gets signaled
        % UE-specific signaling information. The signaled information and
        % where it is located is as follows:
        %   UEsignaling:
        %     - TB_CQI           % CQI used for the transmission of each codeword
        %     - TB_size          % size of the current TB, in bits
        %     - tx_mode          % transmission mode used (SISO, tx diversity, spatial multiplexing)
        %     - rv_idx           % redundancy version index for each codeword
        %   downlink_channel.RB_grid
        %     - user_allocation  % what UE every RB belongs to
        %     - power_allocation % how much power to allocate to each RB,
        %     - n_RB             % RB grid size (frequency)
        %     - sym_per_RB       % number of symbols per RB (12 subcarriers, 0.5ms)
        %     - size_bits        % total size of the RB grid in bits
        %     - numStreams       % maximum number of allowed streams. Resource allocation is described for all of them
        eNodeB_signaling
        
       % Extra tracing options (default options)
       trace_SINR = false;
       
       % average preequalization SNR at current position (averaged over microscopic fading and noise)
       SNR_avg_preequal
       
       % This is an "overall SINR", calculated by summing up all of the
       % signal power and dividing it by the sum of all interfering and
       % noise power. It includes the precoding, so it is not equivalent to
       % the SINR_RS as in the standard. It could nevertheless be used as a
       % similar indicator.
       SINR_RS_not
       
       traffic_model
       lambda = 0;
    end

    methods
        function print(obj)
            if isempty(obj.attached_eNodeB)
                fprintf('User %d, (%d,%d), not attached to an eNodeB\n',obj.id,obj.pos(1),obj.pos(2));
            else
                fprintf('User %d, (%d,%d), eNodeB %d, sector %d\n',obj.id,obj.pos(1),obj.pos(2),obj.attached_eNodeB.id,obj.attached_sector);
            end
            obj.walking_model.print;
        end
        % Move this user according to its settings
        function move(obj)
            new_pos = obj.walking_model.move(obj.pos);
            obj.pos = new_pos;
        end
        % Move this user to where it was the last TTI before according to
        % its settings
        function move_back(obj)
            old_pos = obj.walking_model.move(obj.pos);
            obj.pos = old_pos;
        end
        function UE_in_roi = is_in_roi(a_UE,roi_x_range,roi_y_range)
            % Tells you whether a user in in the Region of Interest (ROI) or not
            % (c) Josep Colom Ikuno, INTHFT, 2008
            % input:    a_UE         ... the UE in question
            %           roi_x_range  ... roi x range. minimum and maximum x coordinates
            %                            which are valid
            %           roi_y_range  ... roi y range. minimum and maximum y coordinates
            %                            which are valid
            % output:   UE_in_roi  ... true or false, whether the UE is inside or not

            UE_pos_x = a_UE.pos(1);
            UE_pos_y = a_UE.pos(2);

            if UE_pos_x<roi_x_range(1) || UE_pos_x>roi_x_range(2)
                UE_in_roi = false;
                return;
            end

            if UE_pos_y<roi_y_range(1) || UE_pos_y>roi_y_range(2)
                UE_in_roi = false;
                return;
            end
            UE_in_roi = true;
        end
        % Starts handover procedures from the currently attached eNodeB to
        % the specified target_eNodeB
        % for now... immediate handover. A proper implementation remains
        % pending.
        function start_handover(obj,target_eNodeB,target_sector)
            % Remove the user from the eNodeB and its scheduler
            obj.attached_eNodeB.sectors(obj.attached_sector).scheduler.remove_UE(obj.id);
            obj.attached_eNodeB.deattachUser(obj);
            % Add the user to the eNodeB and its scheduler
            target_eNodeB.attachUser(obj,target_sector);
            target_eNodeB.sectors(target_sector).scheduler.add_UE(obj.id);
        end

        % Measure whatever needs to be measured and send a feedback to the attached eNodeB
        function send_feedback(obj)
            obj.uplink_channel.send_feedback(obj.feedback);
        end
        
        % Calculates the receiver SINR, which is the metric used to measure link quality
        function link_quality_model(obj,config)
            
            % Get current time
            t = obj.clock.time;
            
            % Get map-dependant parameters for the current user
            interfering_eNodeBs              = obj.attached_eNodeB.neighbors;
            user_macroscopic_pathloss        = obj.downlink_channel.macroscopic_pathloss + obj.penetration_loss + obj.antenna_gain;   % Already includes the TX and RX antenna gain (dB)
            user_macroscopic_pathloss_linear = 10^(0.1*user_macroscopic_pathloss);
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss; % Shadow fading loss (dB)
            user_shadow_fading_loss_linear   = 10^(0.1*user_shadow_fading_loss);

            % Number of codewords, layers, power etc. assigned to this user
            the_RB_grid = obj.downlink_channel.RB_grid;
            tx_mode     = the_RB_grid.tx_mode;

            % Needed so that the dimensions agree
            switch tx_mode
                case 1
                    MIMO = false; % SISO mode
                otherwise
                    MIMO = true;  % MIMO modes
            end
            
            % Get fast fading trace for this subframe
            user_microscale_fading_params      = obj.downlink_channel.fast_fading_pathloss(t,MIMO,tx_mode);
            user_microscale_fading_mode_params = user_microscale_fading_params(tx_mode);

            %% The SINR calculation is done under the following circumstances:
            % Power allocation is done on a per-subframe (1 ms) and RB basis
            % The fast fading trace is given for every 6 subcarriers (every
            % 90 KHz), so as to provide enough samples related to a
            % worst-case-scenario channel length

            % TX power for each layer
            % TX_power_layer_half_RB = the_RB_grid.power_allocation(1,1)/(2*nLayers); % Already divided by 2 (6-subcarrier freq bins) and nLayers NOTE: i'm not sure if it is necessary to divide by nLayers, as already the precoder scales the power
            TX_power_half_RB_data      = reshape([the_RB_grid.power_allocation'; the_RB_grid.power_allocation'],1,[])/(2);     % NOTE: check me! either this line or the one above is correct
            TX_power_half_RB_signaling = reshape([the_RB_grid.power_allocation_signaling'; the_RB_grid.power_allocation_signaling'],1,[])/(2);
            
            TX_power_half_RB = TX_power_half_RB_data + TX_power_half_RB_signaling;
            % TX_power_signaling_half_RB =  TODO: add signaling
            % interference in a nice way!!!!!
            S_dims = size(user_microscale_fading_mode_params.zeta);
            S_dims(2) = 1; % All MATLAB variables have at least 2 dimensions, so not a problem.
            
            % RX power
            if MIMO
                RX_power_half_RB_repmat = repmat(TX_power_half_RB./user_macroscopic_pathloss_linear./user_shadow_fading_loss_linear,S_dims);
                RX_power = RX_power_half_RB_repmat.*user_microscale_fading_mode_params.zeta;
            else
                RX_power = TX_power_half_RB./user_macroscopic_pathloss_linear./user_shadow_fading_loss_linear.*user_microscale_fading_mode_params.zeta.';
            end
            
            % Get interfering eNodeBs
            if ~isempty(interfering_eNodeBs) % no interfering eNodeBs present (single eNodeB simulation)
                interfering_eNodeBs_sectors        = obj.attached_eNodeB.neighbors_eNodeB;
                parent_sites                       = [interfering_eNodeBs_sectors.parent_eNodeB];
                parent_sites_id                    = [parent_sites.id];
                interfering_RB_grids               = [interfering_eNodeBs_sectors.RB_grid];
                interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
                interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
                
                % This line allows us to leave the rest of the code without changes
                interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;
                
                interfering_home_sectors_idx  = 1:length(obj.attached_eNodeB.sectors); % Own sector is later removed by setting its macroscopic pathloss to Inf
                interfering_home_sectors      = obj.attached_eNodeB.sectors(interfering_home_sectors_idx);
                interfering_home_sectors_grid = [interfering_home_sectors.RB_grid];
                interfering_home_sectors_pow_data       = [interfering_home_sectors_grid.power_allocation];
                interfering_home_sectors_pow_signaling  = [interfering_home_sectors_grid.power_allocation_signaling];
                
                % This line allows us to leave the rest of the code without changes
                interfering_home_sectors_pow = interfering_home_sectors_pow_data + interfering_home_sectors_pow_signaling;
                
                interfering_home_sectors_id   = [interfering_home_sectors.eNodeB_id];
                interfering_sectors_id        = [interfering_eNodeBs_sectors.eNodeB_id];
                
                % Add the neighboring sectors
                interfering_sectors_id_all    = [interfering_home_sectors_id interfering_sectors_id];
                
                if ~MIMO
                    microscale_interfering_thetas = user_microscale_fading_mode_params.theta(:,interfering_sectors_id_all).';
                else
                    microscale_interfering_thetas = user_microscale_fading_mode_params.theta(:,:,interfering_sectors_id_all,:);
                end
                SINR_interf_dims = size(microscale_interfering_thetas);
                SINR_interf_dims_repmat = SINR_interf_dims;
                SINR_interf_dims_repmat(2) = 1;
                if length(SINR_interf_dims_repmat) > 2
                    SINR_interf_dims_repmat(3) = 1;
                end
                
                % Get assigned interfering power
                if config.feedback_channel_delay~=0
                    % Take scheduled power
                    interf_power_all    = [interfering_home_sectors_pow interfering_power_allocations];
                    interf_power_all_RB = kron(interf_power_all,[1;1])/2';
                else
                    % Turn on all interferers
                    interf_power_all_RB = repmat([[interfering_eNodeBs_sectors.max_power] [interfering_home_sectors.max_power]]/SINR_interf_dims(2),[SINR_interf_dims(2) 1]);
                end
                if MIMO
                    interf_power_all_RB_repmat = reshape(repmat(interf_power_all_RB,SINR_interf_dims_repmat),SINR_interf_dims); % Also valid for the case where more than one rank is used
                else
                    interf_power_all_RB_repmat = interf_power_all_RB.';
                end
                
                % Macro and Shadow fading pathloss. Site and eNodeB dependent
                interferingEnodeBids = [obj.attached_eNodeB.id interfering_eNodeBs.id]; % Site list
                home_sectors = [obj.attached_eNodeB.sectors];
                home_parent  = [home_sectors.parent_eNodeB];
                interferingSiteIds = [[home_parent.id] parent_sites_id];
                interfering_macroscopic_pathloss_eNodeB = obj.downlink_channel.interfering_macroscopic_pathloss(interferingEnodeBids) + obj.penetration_loss + obj.antenna_gain;
                
                % Take out the own sector from the interference power matrix by setting the macroscopic pathloss to Inf
                interfering_macroscopic_pathloss_eNodeB(obj.attached_sector) = Inf; % The list begins with the home sectors
                
                % Continue with macro and shadow fading data retrieving
                interfering_macroscopic_pathloss_eNodeB_linear = 10.^(0.1*interfering_macroscopic_pathloss_eNodeB);
                [interfering_shadow_fading_loss shadow_fading_is_set] = obj.downlink_channel.interfering_shadow_fading_pathloss(interferingSiteIds);
                
                % For the case that no shadow fading  considered (using network planning tool)
                if ~shadow_fading_is_set
                    interfering_shadow_fading_loss = repmat(interfering_shadow_fading_loss,[size(interfering_macroscopic_pathloss_eNodeB_linear,1) 1]);
                else
                    % Do nothing, as the vetor already has the correct length
                end
                interfering_shadow_fading_loss_linear = 10.^(0.1*interfering_shadow_fading_loss);
                
                temp_macro_mat = interfering_macroscopic_pathloss_eNodeB_linear';
                temp_macro_mat = temp_macro_mat(ones(SINR_interf_dims(2),1),:);
                temp_shadow_mat = interfering_shadow_fading_loss_linear';
                temp_shadow_mat = temp_shadow_mat(ones(SINR_interf_dims(2),1),:);
                if MIMO
                    interf_macro_fading  = reshape(repmat(temp_macro_mat,SINR_interf_dims_repmat),SINR_interf_dims);
                    interf_shadow_fading = reshape(repmat(temp_shadow_mat,SINR_interf_dims_repmat),SINR_interf_dims);
                else
                    interf_macro_fading  = temp_macro_mat.';
                    interf_shadow_fading = temp_shadow_mat.';
                end
            else
               interf_macro_fading = 1;
               interf_shadow_fading = 1;
               microscale_interfering_thetas = 0; % set interference power equal to zero
            end
            % Calculate thermal noise
            thermal_noise_watts_per_RB = 10^(0.1*(obj.downlink_channel.thermal_noise_dBW_RB + obj.receiver_noise_figure));
            
            % Calculate average preequalization SNR
            % This is a total SNR, the same as in the Link Level Simulator (the channel is normalized to nTX, so this is also necessary here!)
            obj.SNR_avg_preequal = 10*log10(TX_power_half_RB(1)./(user_macroscopic_pathloss_linear*user_shadow_fading_loss_linear*thermal_noise_watts_per_RB/2));%*config.nTX));

            if ~MIMO
                % SINR calculation
                noise_plus_inter_layer_power = user_microscale_fading_mode_params.psi.*thermal_noise_watts_per_RB/2;
                
                if ~isempty(interfering_eNodeBs)
                    % Also works for more than one rank (i.e. extra dimension)
                    interfering_rx_power = squeeze(sum(interf_power_all_RB_repmat./interf_macro_fading./interf_shadow_fading.*microscale_interfering_thetas,1));
                    Interference_plus_noise_power = noise_plus_inter_layer_power + interfering_rx_power.';
                else
                    Interference_plus_noise_power = noise_plus_inter_layer_power;
                end
                SINR_linear = RX_power ./ Interference_plus_noise_power.'; % Divide thermal noise by 2: Half-RB frequency bins
            else
                noise_plus_inter_layer_power = user_microscale_fading_mode_params.chi.*RX_power + user_microscale_fading_mode_params.psi.*thermal_noise_watts_per_RB/2; % Divide thermal noise by 2: Half-RB frequency bins
                if ~isempty(interfering_eNodeBs)
                    % Also works for more than one rank (i.e. extra dimension)
                    interfering_rx_power = squeeze(sum(interf_power_all_RB_repmat./interf_macro_fading./interf_shadow_fading.*microscale_interfering_thetas,3));
                    Interference_plus_noise_power = noise_plus_inter_layer_power + interfering_rx_power;
                else
                    Interference_plus_noise_power = noise_plus_inter_layer_power;
                end
                SINR_linear = RX_power./Interference_plus_noise_power;
            end
            
            % Calculation of "overal SINR". Which should be similar to
            % the RS SINR (although differente, as the RS SINR should not
            % include precoding and this does)
            RX_power_all                      = squeeze(sum(sum(RX_power,2),1));
            Interference_plus_noise_power_all = squeeze(sum(sum(Interference_plus_noise_power,2),1));
            obj.SINR_RS_not                   = 10*log10(RX_power_all(1)/Interference_plus_noise_power_all(1)); % since it includes the precoding, I cannot use it as RS SINR
            
            % Calculation of the post-equalization symbols SINR
            SINR_dB = 10*log10(SINR_linear);
            
            % Calculate and save feedback, as well as the measured SINRs
            obj.calculate_feedback(config,tx_mode,SINR_linear,SINR_dB);
        end
        
        % Calculate the feedback values based on the input. This function
        % is called from the link quality model and is separated for
        % convenience and readability. The results of the feedback
        % calculation are stored in the following variables:
        % - obj.feedback.CQI:              CQI feedback
        % - obj.feedback.RI:               Rank Indicator feedback (when applicable)
        % - obj.link_quality_model_output: SINR values
        %
        % As input parameters you have one SINR per RB
        % (SINRs_to_map_to_CQI) or all of the SINRs the SL simulator
        % traces, which are currently two per RB (SINR_dB)
        function calculate_feedback(obj,config,tx_mode,SINR_linear,SINR_dB)
            % Take a subset of the SINRs for feedback calculation
            % For SM we send 2 CQIs, one for each of the codewords (which in the 2x2
            % case are also the layers). For TxD, both layers have the same SINR
            % The CQI is calculated as a linear averaging of the SINRs in
            % dB. This is done because like this the Tx has an "overall
            % idea" of the state of the RB, not just a sample of it.
            switch tx_mode
                case 1 % SISO
                    SINRs_to_map_to_CQI = (SINR_dB(1:2:end)+SINR_dB(2:2:end))/2;
                    obj.link_quality_model_output.SINR_dB     = SINR_dB;
                    obj.link_quality_model_output.SINR_linear = SINR_linear;
                case 2 % TxD
                    % Both layers have the same SINR
                    SINRs_to_map_to_CQI = (SINR_dB(1,1:2:end)+SINR_dB(1,2:2:end))/2;
                    obj.link_quality_model_output.SINR_dB     = SINR_dB(1,:);
                    obj.link_quality_model_output.SINR_linear = SINR_linear(1,:);
                case {3,4} % OLSM, CLSM
                    SINRs_to_map_to_CQI = (SINR_dB(:,1:2:end,:)+SINR_dB(:,2:2:end,:))/2;
                    obj.link_quality_model_output.SINR_dB     = SINR_dB;
                    obj.link_quality_model_output.SINR_linear = SINR_linear;
                otherwise
                    error('TX mode not yet supported');
            end
            
            % Send as feedback the CQI for each RB.
            % Flooring the CQI provides much better results than
            % rounding it, as by rounding it to a higher CQI you will
            % very easily jump the BLER to 1. The other way around it
            % will jump to 0.
            if obj.unquantized_CQI_feedback
                CQIs = obj.CQI_mapper.SINR_to_CQI(SINRs_to_map_to_CQI);
            else
                CQIs = floor(obj.CQI_mapper.SINR_to_CQI(SINRs_to_map_to_CQI));
            end
            if (tx_mode==3) || (tx_mode==4) % Rank decision for SM
                % Decide based on the number of transmitted data bits for a rank value
                bits2 = zeros(size(CQIs,3),1);
                for RI_idx = 1:size(CQIs,3)
                    % Layer-to-codeword mapping (for now just up to 2 layers: easy)
                    SINRs_current_cw = SINRs_to_map_to_CQI(1:RI_idx,:,RI_idx);
                    for CW_idx = 1:RI_idx % In the future this will need to be fixed to be able to accomodate other layer mappings. but for now it is OK
                        [SINR_av_lin SINR_av_dB] = obj.SINR_averager.average(SINRs_current_cw(CW_idx,:),1:15,true); % Input directly in dB
                        CQI_temp = floor(obj.CQI_mapper.SINR_to_CQI(SINR_av_dB));
                        CQI_layer = max(CQI_temp(~(CQI_temp - (1:15)')));
                        if CQI_layer
                            bits2(RI_idx) = bits2(RI_idx) + 8*round(1/8*config.CQI_params(CQI_layer).modulation_order * config.CQI_params(CQI_layer).coding_rate_x_1024/1024 * config.sym_per_RB_nosync * config.N_RB*2)-24;
                        end
                    end
                end
                [C,rank_i]       = max(bits2);              % Choose the RI for which the number of bits is maximized
                obj.feedback.CQI = CQIs(1:rank_i,:,rank_i); % Feedback of the chosen rank
                obj.feedback.RI  = rank_i;
            else
                obj.feedback.CQI = CQIs;
                obj.feedback.RI  = 1;
            end
            obj.feedback.tx_mode = tx_mode;
        end
        
        % Evaluate whether this TB arrived correctly by using the data from
        % the link quality model and feeding it to the link performance
        % model (BLER curves)
        function link_performance_model(obj)
            
            % Get RB grid
            % the_RB_grid   = obj.RB_grid;
            
            % Get SINRs from the link quality model. Only the dB (not
            % linear) are needed.
            SINR_dB       = obj.link_quality_model_output.SINR_dB;
                      
            % Calculate TB SINR
            TB_CQI       = obj.eNodeB_signaling.TB_CQI;
            user_RBs     = obj.eNodeB_signaling.assigned_RB_map;
            assigned_RBs = obj.eNodeB_signaling.num_assigned_RBs;
            tx_mode      = obj.eNodeB_signaling.tx_mode;
            nLayers      = obj.eNodeB_signaling.nLayers;
            nCodewords   = obj.eNodeB_signaling.nCodewords;
            rv_idxs      = obj.eNodeB_signaling.rv_idx;
            TB_size      = obj.eNodeB_signaling.TB_size;
            
            % Preallocate variables to store in trace
            TB_SINR_dB = zeros(1,nCodewords);
            BLER       = zeros(1,nCodewords);
            
            % Needed for non-full-buffer simulations
            N_used_bits  = obj.eNodeB_signaling.N_used_bits;
            packet_parts = obj.eNodeB_signaling.packet_parts;
            
            % NOTE: This needs to be changed into a per-layer SINR. A layer-to-stream translation will be needed
            UE_TB_SINR_idxs = logical(kron(user_RBs',[1 1]));
            % Setting this Kronecker product to [1 0] would mean that you only take one SC per RB and repeat it. i.e. not use all of
            % your SCs in the trace.

            % Set feedback for all streams
            if assigned_RBs~=0
                % Not all of the dimensions are needed
                switch tx_mode
                    case {1,2}
                        % SIXO, TxD
                    case {3,4}
                        % OLSM, CLSM
                        
                        % Take only the SINRs of the layers that are needed
                        SINR_dB = SINR_dB(1:nLayers,:,nLayers);
                        
                        % Layer mapping
                        % NOTHING (for now) -> 2 TX antennas does not yet need a proper mapping
                        % Convert the shape of the SINR_dB vector, as well
                        % as the mapping
                    otherwise
                        error('Mode not supported');
                end
                
                for l_=1:nLayers
                    layer_SINRs       = SINR_dB(l_,:);
                    UE_TB_SINRs_layer = layer_SINRs(UE_TB_SINR_idxs);
                    [TB_SINR_lin TB_SINR_dB(l_)]    = obj.SINR_averager.average(UE_TB_SINRs_layer,TB_CQI(l_),true);
                    BLER(l_)          = obj.BLER_curves.get_BLER(TB_CQI(l_),TB_SINR_dB(l_));
                end
                % Receive
                ACK     = BLER<rand(1,nCodewords);
            else
                % Dummy results
                TB_SINR_dB = [];
                ACK        = false(1,nCodewords);
            end
            
            % Needed for non-full-buffer simulations
            if ~obj.traffic_model.is_fullbuffer
                for cw_ = 1:nCodewords
                    if ACK(cw_)
                        if strcmp(obj.traffic_model.type,'voip') || strcmp(obj.traffic_model.type,'video') || strcmp(obj.traffic_model.type,'gaming')
                            for pp = 1:length(packet_parts{cw_}) % acknowledge all packet parts and remove them from the buffer
                                if packet_parts{cw_}(pp).data_packet_id
                                    packet_ind = obj.traffic_model.get_packet_ids == packet_parts{cw_}(pp).data_packet_id;
                                    [packet_done,packet_id] = obj.traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(packet_parts{cw_}(pp).id,false);
                                    if packet_done && packet_id
                                        obj.traffic_model.remove_packet(packet_id,true);
                                    end
                                end
                            end
                        end
                    else
                        N_used_bits(cw_) = 0;
                        if strcmp(obj.traffic_model.type,'voip') || strcmp(obj.traffic_model.type,'video') || strcmp(obj.traffic_model.type,'gaming')
                            for pp = 1:length(packet_parts{cw_})
                                if packet_parts{cw_}(pp).data_packet_id
                                    packet_ind = obj.traffic_model.get_packet_ids == packet_parts{cw_}(pp).data_packet_id;
                                    [packet_done,packet_id] = obj.traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(packet_parts{cw_}(pp).id,false);
                                    if packet_done && packet_id
                                        obj.traffic_model.remove_packet(packet_id,false);
                                    end
                                end
                            end
                        end
                    end
                end
            end
                      
            % Add BLER/ACK feedback to the CQI and RI feedback
            if assigned_RBs~=0
                obj.feedback.UE_scheduled = true;
                obj.feedback.nCodewords   = nCodewords;
                obj.feedback.TB_size      = TB_size;
                obj.feedback.BLER         = BLER;
                obj.feedback.ACK          = ACK;
            else
                obj.feedback.UE_scheduled = false;
                obj.feedback.nCodewords   = 0;
                obj.feedback.TB_size      = 0;
                obj.feedback.BLER         = NaN;
                obj.feedback.ACK          = false;
            end
            
            % Optional traces
            if obj.trace_SINR
                extra_traces{1} = SINR_dB;
                extra_traces{2} = obj.SNR_avg_preequal;
            else
                extra_traces{1} = [];
                extra_traces{2} = [];
            end
            
            % Store trace of the relevant information
            tti_idx = obj.clock.current_TTI;
            
            % Store trace
            obj.trace.store(...
                nCodewords,...
                obj.feedback.CQI,...
                obj.attached_eNodeB.id,...
                obj.attached_sector,...
                obj.pos,tti_idx,...
                assigned_RBs,...
                ACK,...
                TB_CQI,...
                TB_size,...
                BLER,...
                TB_SINR_dB,...
                N_used_bits,...
                obj.SINR_RS_not,...
                extra_traces);
        end
    end
end
