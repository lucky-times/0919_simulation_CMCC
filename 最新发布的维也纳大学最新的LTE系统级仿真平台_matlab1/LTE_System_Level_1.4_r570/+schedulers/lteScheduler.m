classdef lteScheduler < handle
    % Implements common methods needed by any implementation of an LTE
    % scheduler (eg. Round Robin, Best CQI...)
    % (c) Josep Colom Ikuno, INTHFT, 2009
    
    properties
        % Scheduler name
        name = 'Generic LTE scheduler superclass';
        
        % Where this scheduler is attached
        attached_eNodeB_sector
        % Copy of the CQI tables
        CQI_range
        CQI_tables
        CQIs_efficiency % In order to avoid getting this every TTI
        % The algorithm used to average the several SINRs into a single TB SINR
        SINR_averager
        % BLER data and CQI mapping
        BLER_curves
        CQI_mapper
        % To control power allocation
        max_power
        
        % Target system BLER
        target_BLER = 0.1;
        
        % Set this option to false if you are doing SL-LL validation: if not, this leads to problems with the simulated
        % BLER. For CQI 1 one could not simulate BLERs worse than 0.1 as those RBs are just skipped
        skip_null_CQIs = true;
        
        % Where trace data is saved
        trace
        
        % Genie information (eg: all of the eNodeBs and UEs)
        genie
        
        % UE trace: from where the historical information (ie. UE throughput) is extracted
        UE_traces
        % In order to know how many TTIs to skip when looking for throughput
        % (in reality you would not know instantly whether the TB was received correctly or not)
        feedback_delay_TTIs
        
        % Network clock. Tells the scheduler in what TTI he is
        clock
        
        % Other things
        av_const
        fairness
        k
        d
        MUMIMO
        overhead_ref
        overhead_sync
    end
    
    methods(Abstract)
        schedule_users(obj,RB_grid,attached_UEs,UE_feedback)
        add_UE(obj,UE_id)
        remove_UE(obj,UE_id)
    end
    
    methods
        
        % Class constructor
        function obj = lteScheduler(scheduler_params,attached_eNodeB_sector)
            obj.attached_eNodeB_sector = attached_eNodeB_sector;
            CQI_range = LTE_common_get_CQI_params('range');
            obj.CQI_tables = LTE_common_get_CQI_params(CQI_range(1)); % To initialize the struct
            for i_=CQI_range(1):CQI_range(2)
                obj.CQI_tables(i_) = LTE_common_get_CQI_params(i_);
            end
            obj.CQI_range = CQI_range(1):CQI_range(2);
            obj.CQIs_efficiency = [0 [obj.CQI_tables.efficiency]]; % Take note of CQI 0 also
            obj.max_power = scheduler_params.max_power;
            obj.clock = attached_eNodeB_sector.parent_eNodeB.clock;
            obj.av_const = scheduler_params.av_window;
            obj.fairness = scheduler_params.fairness;
            obj.k = scheduler_params.k;
            obj.d = scheduler_params.d;
            obj.overhead_ref = scheduler_params.overhead_ref;
            obj.overhead_sync = scheduler_params.overhead_sync;
        end
        
        % Print some info
        function print(obj)
            fprintf('%s\n',obj.name);
        end
        
        % Find the optimum CQI values for a set of N codewords. It is
        % assumed that the CQI vector/matrix is NOT empty
        function [assigned_CQI predicted_BLER predicted_SINR] = get_optimum_CQIs(obj,CQIs_to_average_all,nCodewords)
            % Preallocation
            assigned_CQI   = zeros(nCodewords,1);
            predicted_BLER = zeros(nCodewords,1);
            predicted_SINR = zeros(nCodewords,1);
            
            UE_estimated_SINR_dB_all = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average_all);
            for cw_=1:nCodewords
                % don't use the SINR value on the lower bound of a CQI interval
                % but something inbetween (otherwise this is too conservative)
                % - as in the Link Level Simulator (heuristically optimized)
                %             SINR_temp1 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average);
                %             SINR_temp2 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average+1);
                % %
                %             UE_estimated_SINR_dB = SINR_temp1+(SINR_temp2-SINR_temp1)/2.5;
                %             UE_estimated_SINR_linear = 10.^(0.1*UE_estimated_SINR_dB);
                [averaged_SINR_MCS_dependent_lin averaged_SINR_MCS_dependent_dB] = obj.SINR_averager.average(UE_estimated_SINR_dB_all(:,cw_),obj.CQI_range,true);
                predicted_BLERs                = obj.BLER_curves.get_BLER_CQI(obj.CQI_range,averaged_SINR_MCS_dependent_dB);
                
                % Objective is the closest smaller or equal to 10% BLER (BLER 0 is preferred to BLER 1)
                if predicted_BLERs(end) == 0
                    % Case of a very good channel
                    cw_assigned_CQI = 15;
                elseif predicted_BLERs(1) >= obj.target_BLER
                    % Case of a bad channel
                    cw_assigned_CQI = 1;
                else
                    % Case in the middle
                    abs_diffs = predicted_BLERs-obj.target_BLER;
                    abs_diffs = round(abs_diffs*1000)/1000; % To avoid small statistical mistakes in the BLER plots. No change assuming that the target BLER is in the order of 10%
                    cw_assigned_CQI = find(abs_diffs<=0,1,'last');
                end
                assigned_CQI(cw_)   = cw_assigned_CQI;
                predicted_BLER(cw_) = predicted_BLERs(cw_assigned_CQI);
                predicted_SINR(cw_) = averaged_SINR_MCS_dependent_dB(cw_assigned_CQI);     % effective logarithmic SINR
            end
        end
        
        % For a set of SINR values, this functions averages them using the
        % scheduler's SINR averaging algorithm and outputs the MCS that
        % most closely approaches the target BLER (set at 0.1)
        function [assigned_CQI predicted_BLER effective_SINR] = get_optimum_CQI(obj,UE_estimated_SINR_dB)
            % don't use the SINR value on the lower bound of a CQI interval
            % but something inbetween (otherwise this is too conservative)
            % - as in the Link Level Simulator (heuristically optimized)
            %             SINR_temp1 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average);
            %             SINR_temp2 = obj.CQI_mapper.CQI_to_SINR(CQIs_to_average+1);
            % %
            %             UE_estimated_SINR_dB = SINR_temp1+(SINR_temp2-SINR_temp1)/2.5;
            %             UE_estimated_SINR_linear = 10.^(0.1*UE_estimated_SINR_dB);
            
            [averaged_SINR_MCS_dependent_lin averaged_SINR_MCS_dependent_dB] = obj.SINR_averager.average(UE_estimated_SINR_dB,obj.CQI_range,true);
            predicted_BLERs                = obj.BLER_curves.get_BLER_CQI(obj.CQI_range,averaged_SINR_MCS_dependent_dB);
            
            % Objective is the closest smaller or equal to 10% BLER (BLER 0 is preferred to BLER 1)
            if predicted_BLERs(end) == 0
                % Case of a very good channel
                assigned_CQI = 15;
            elseif predicted_BLERs(1) >= obj.target_BLER
                % Case of a bad channel
                assigned_CQI = 1;
            else
                % Case in the middle
                abs_diffs = predicted_BLERs-obj.target_BLER;
                abs_diffs = round(abs_diffs*1000)/1000; % To avoid small statistical mistakes in the BLER plots. No change assuming that the target BLER is in the order of 10%
                assigned_CQI = find(abs_diffs<=0,1,'last');
            end
            predicted_BLER = predicted_BLERs(assigned_CQI);
            effective_SINR = averaged_SINR_MCS_dependent_dB(assigned_CQI);     % effective logarithmic SINR
        end
        
        function CQI_spectral_efficiency = get_spectral_efficiency(obj,CQI_matrix)
            % Returns the spectral efficiencies related to the CQIs in the
            % matrix. In case of unquantized feedback, it floors the CQI
            % values to the nearest usable CQI
            CQI_idx                 = floor(CQI_matrix)+1; % Get CQI indexes
            CQI_idx(CQI_idx<1)      = 1;       % Clip extremes
            CQI_idx(CQI_idx>16)     = 16;
            sizes_CQI_feedback      = size(CQI_idx);
            CQI_idx_vector          = uint16(CQI_idx(:)); % To avoid an error caused by Matlab sometimes saying those were not integers
            CQI_spectral_efficiency = reshape(obj.CQIs_efficiency(CQI_idx_vector),sizes_CQI_feedback);
        end
        
        function UE_resource_assignment = get_max_UEs(obj,metric_matrix,UE_assignment)
            % Scans a nUExM matrix and returns for each column the index of
            % the highest value. In case more than one UE has the same
            % metric value, one of these ones is randomly selected. Assume
            % all values to be non-negative and 0 as "out-of-range"
            max_metric = max(metric_matrix,[],1);
            UE_resource_assignment = zeros(length(max_metric),1);
            for rb_=1:length(max_metric)
                if max_metric(rb_)~=0
                    candidates = find(metric_matrix(:,rb_)==max_metric(rb_));
                    index = ceil(rand*length(candidates));
                    if index==0
                        index = 1;
                    end
                    UE_resource_assignment(rb_) = UE_assignment(candidates(index)); % Choose a random UE from the set
                end
            end
        end
        
        function [N_assigned_RBs CQIs_to_average_all UE_scheduled new_UE_RB_map] = filter_out_zero_RBs_and_get_CQIs(obj,RB_grid,nCodewords,UE_CQI_feedback,current_UE)
            % Do not use RBs with a CQI of 0 (they are lost)
            UE_RBs = RB_grid.user_allocation==current_UE.id;
            if obj.skip_null_CQIs
                if nCodewords == 1
                    zero_CQIs     = (UE_CQI_feedback(:,1)<1);
                    non_valid_RBs = UE_RBs & zero_CQIs;  % RBs to filter out
                    valid_RBs     = UE_RBs & ~zero_CQIs; % CQIs that will be averaged
                    RB_grid.user_allocation(non_valid_RBs) = 0;
                    CQIs_to_average_all = UE_CQI_feedback(valid_RBs);
                else
                    % For the case where more than 1 codewords are sent, all CWs must have a CQI >0
                    zero_CQIs     = sum(UE_CQI_feedback<1,2)>=1;
                    non_valid_RBs = UE_RBs & zero_CQIs; % RBs to filter out
                    valid_RBs     = UE_RBs & ~zero_CQIs; % CQIs that will be averaged
                    RB_grid.user_allocation(non_valid_RBs) = 0;
                    CQIs_to_average_all = UE_CQI_feedback(valid_RBs,:);
                end
                new_UE_RB_map  = valid_RBs;
                N_assigned_RBs = sum(valid_RBs);
            else
                if nCodewords == 1
                    CQIs_to_average_all = UE_CQI_feedback(UE_RBs);
                else
                    CQIs_to_average_all = UE_CQI_feedback(UE_RBs,:);
                end
                new_UE_RB_map  = UE_RBs;
                N_assigned_RBs = sum(UE_RBs);
            end
            
            if isempty(CQIs_to_average_all)
                UE_scheduled   = false;
            else
                UE_scheduled = true;
            end
        end
        
        function schedule_users_common(obj,RB_grid,attached_UEs,UE_feedback,current_TTI,tx_mode)
            % NOTE: since for the zero-delay case, RI and nLayers of the
            % transmission may not match, nCodewords is not used here
            
            % Common tasks for all schedulers
            RB_grid_size_bits  = 0;
            max_codewords      = 2;
            predicted_UE_BLERs = NaN(max_codewords,length(attached_UEs));
            assigned_UE_RBs    = zeros(1,length(attached_UEs));
            
            % Homogeneous power allocation
            if ~isempty(attached_UEs)
                RB_grid.power_allocation(:) = obj.max_power / RB_grid.n_RB;
            end
            
            % Continue UE common scheduling procedures
            for u_=1:length(attached_UEs)
                current_UE = attached_UEs(u_);
                
                % If the feedback is present: process normally
                if UE_feedback.feedback_received(u_)
                    
                    tx_mode = UE_feedback.tx_mode(u_);
                    if (tx_mode == 4) || (tx_mode == 3)  % in SM modes use the RI value also
                        nLayers = UE_feedback.RI(u_);
                        nCodewords = min(2,nLayers); % To be changed for the 4+ layer cases
                    else
                        nLayers    = 1;
                        nCodewords = 1;
                    end
                    
                    UE_CQI_feedback = UE_feedback.CQI(:,1:nLayers,u_);
                    
                    % Do not use RBs with a CQI of 0 (they are lost).
                    % This function also averages the CQIs and assigns an overall TB CQI with predicted BLER < 10%
                    [num_assigned_RB CQIs_to_average_all UE_scheduled RB_map] = obj.filter_out_zero_RBs_and_get_CQIs(RB_grid,nCodewords,UE_CQI_feedback,current_UE);
                    
                    if UE_scheduled
                        % Simplified this piece of code by using the superclass, as all types of scheduler will to make use of it.
                        [assigned_CQI predicted_UE_BLERs(1:nCodewords,u_) estimated_TB_SINR] = obj.get_optimum_CQIs(CQIs_to_average_all,nCodewords);
                        % Signal down the user CQI assignment
                        attached_UEs(u_).eNodeB_signaling.assigned_RB_map = RB_map;
                        attached_UEs(u_).eNodeB_signaling.tx_mode         = tx_mode;
                        attached_UEs(u_).eNodeB_signaling.TB_CQI          = assigned_CQI;
                        attached_UEs(u_).eNodeB_signaling.nCodewords      = nCodewords;
                        attached_UEs(u_).eNodeB_signaling.nLayers         = nLayers;
                        attached_UEs(u_).eNodeB_signaling.genie_TB_SINR   = estimated_TB_SINR;
                    else
                        attached_UEs(u_).eNodeB_signaling.assigned_RB_map = [];
                        attached_UEs(u_).eNodeB_signaling.tx_mode         = 0;
                        attached_UEs(u_).eNodeB_signaling.TB_CQI          = 0;
                        attached_UEs(u_).eNodeB_signaling.nCodewords      = 0;
                        attached_UEs(u_).eNodeB_signaling.nLayers         = 0;
                        attached_UEs(u_).eNodeB_signaling.genie_TB_SINR   = NaN;
                    end
                    
                else
                    % If the feedback is not present: assign a default CQI value of 1 with rank one.
                    UE_scheduled    = true;
                    RB_map          = RB_grid.user_allocation==current_UE.id;
                    num_assigned_RB = sum(RB_map);
                    nCodewords      = 1;
                    nLayers         = 1;
                    % Signal down the user CQI assignment
                    attached_UEs(u_).eNodeB_signaling.assigned_RB_map      = RB_map;
                    attached_UEs(u_).eNodeB_signaling.tx_mode              = tx_mode;
                    attached_UEs(u_).eNodeB_signaling.TB_CQI(1:nCodewords) = 1;
                    attached_UEs(u_).eNodeB_signaling.nCodewords           = nCodewords;
                    attached_UEs(u_).eNodeB_signaling.nLayers              = nLayers;
                    attached_UEs(u_).eNodeB_signaling.genie_TB_SINR        = NaN;
                    predicted_UE_BLERs(u_) = 0; % Dummy value to avoid a NaN
                end
                
                if UE_scheduled
                    TB_CQI_params    = obj.CQI_tables(attached_UEs(u_).eNodeB_signaling.TB_CQI);
                    modulation_order = [TB_CQI_params.modulation_order];
                    coding_rate      = [TB_CQI_params.coding_rate_x_1024]/1024;
                    
                    if mod(current_TTI-1,5) % TB without sync symbols
                        % The factor of two is because there are two time slots per subframe, 24 CRC bits are attached
                        % Segmentation prior to channel coding NOT taken into account.
                        TB_size_bits = max(8*round(1/8*(RB_grid.sym_per_RB_nosync .* num_assigned_RB .* modulation_order .* coding_rate * 2))-24,0);
                    else % TB with sync symbols
                        sync_pos = false(size(RB_grid.user_allocation));
                        sync_pos(floor(length(sync_pos)/2)-2:floor(length(sync_pos)/2)+3) = true;
                        sync_RBs = sum((RB_grid.user_allocation==attached_UEs(u_).id) .* sync_pos);
                        non_sync_RBs = sum(RB_grid.user_allocation==attached_UEs(u_).id)-sync_RBs;
                        TB_size_bits = max(8*round(1/8*(RB_grid.sym_per_RB_sync .* sync_RBs + RB_grid.sym_per_RB_nosync * non_sync_RBs) .* modulation_order .* coding_rate * 2)-24,0);
                    end
                else
                    num_assigned_RB = 0;
                    TB_size_bits = 0;
                end
                attached_UEs(u_).eNodeB_signaling.num_assigned_RBs = num_assigned_RB;
                attached_UEs(u_).eNodeB_signaling.TB_size = TB_size_bits;
                attached_UEs(u_).eNodeB_signaling.rv_idx = 0;
                RB_grid_size_bits = RB_grid_size_bits + TB_size_bits;
                
                assigned_UE_RBs(u_) = num_assigned_RB;
                
                for cw_ = 1:length(TB_size_bits)
                    if TB_size_bits(cw_) ~= 0
                        if (~attached_UEs(u_).traffic_model.is_fullbuffer) && (strcmp(attached_UEs(u_).traffic_model.type,'voip') || strcmp(attached_UEs(u_).traffic_model.type,'video') || strcmp(attached_UEs(u_).traffic_model.type,'gaming'))
                            packet_parts = attached_UEs(u_).traffic_model.decrease_packets(TB_size_bits(cw_));
                            if ~isempty(packet_parts)
                                attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = sum(packet_parts.get_size);
                                %                                     obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts = packet_parts;
                            else
                                attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = 0;
                            end
                            attached_UEs(u_).eNodeB_signaling.packet_parts{cw_} = packet_parts;
                        else
                            attached_UEs(u_).traffic_model.decrease_packets(TB_size_bits(cw_));
                            attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) = min(attached_UEs(u_).traffic_model.get_buffer_length,TB_size_bits(cw_));
                        end
                    else
                        attached_UEs(u_).eNodeB_signaling.N_used_bits(cw_) =  0;
                        attached_UEs(u_).eNodeB_signaling.packet_parts{cw_} = [];
                    end
                end
            end
            
            RB_grid.size_bits = RB_grid_size_bits;
            
            % TODO: HARQ handling, #streams decision and tx_mode decision.
            
            % Store trace
            TTI_idx = obj.attached_eNodeB_sector.parent_eNodeB.clock.current_TTI;
            obj.trace.store(TTI_idx,mean(assigned_UE_RBs),mean(predicted_UE_BLERs(isfinite(predicted_UE_BLERs))));
        end
        
        function TP = compute_av_throughput(obj,u_,UE_feedback,TTI_to_read)
            UE_id = UE_feedback.UE_id(u_);
            if UE_id
                TP = sum(obj.UE_traces(UE_id).avg_throughput(:,TTI_to_read))*10^-3; % Mean throughput, averaged with an exponential window
            else
                TP = 0;
            end
        end
        
        function [c,user_ind] = get_efficiency(obj,N_UE,N_RB,UE_feedback)
            a = zeros(N_UE,1);
            b = zeros(N_UE,1);
            c = zeros(N_UE,N_RB);
            
            for u_ = 1:N_UE
                CQI_bar = max(UE_feedback.CQI(u_,:))+1;
                a(u_) = obj.k(CQI_bar);
                b(u_) = obj.d(CQI_bar);
            end
            
            user_ind = randperm(N_UE);
            for rb = 1:N_RB
                c_count = 0;
                for uu = user_ind
                    c_count = c_count+1;
                    c(c_count,rb) = a(uu)* UE_feedback.CQI(uu,rb)+b(uu);
                end
            end
        end
    end
    
end

