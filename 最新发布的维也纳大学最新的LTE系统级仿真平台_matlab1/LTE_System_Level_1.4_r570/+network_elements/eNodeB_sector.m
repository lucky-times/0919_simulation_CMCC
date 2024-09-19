classdef eNodeB_sector < handle
    % Defines an eNodeB's sector
    % (c) Josep Colom Ikuno, INTHFT, 2008

    properties
        % Sector id inside the site
        id
        % Sector id for the whole sector (eNodeB) set
        eNodeB_id
        % eNodeB to which this sector belongs
        parent_eNodeB
        % Sector antenna's azimuth
        azimuth
        % Sector antenna
        antenna
        % Users attached to this eNodeB. Array of UEs stored as a
        % linked list. This points to the first user
        users = [];
        % number of attached UEs
        attached_UEs = 0;
        % eNodeB sector maximum transmit power for data, in Watts
        max_power
        % Power dedicated to signaling. Counted as always in use
        signaling_power
        % This sector's scheduler
        scheduler
        % This sector's currently used resource block assignment grid
        RB_grid
        % Number of antennas
        nTX
        
        % Controls whether the eNodeB is always radiating power (dafault and worse-case scenario) or no power is used when no UEs are attached
        always_on = true;
        
        % the last received feedback and a list of attached UEs
        last_received_feedback
        UEs_attached_last_TTI
        
        % trace that stores the received feedbacks
        feedback_trace
        % trace that stores the RB assignments
        sector_trace
        
        % Configuration option for zero-delay CQI feedback
        zero_delay_feedback = false;
        
        % this parameter allows you to send unquantized feedback. It can
        % serve to assert how good a certain CQI mapping approaches the
        % target BLER in comparison to directly knowing the SINR
        unquantized_CQI_feedback = false;
        
        % extra params that may not be always used (mainly extra information from the network planning tool)
        transmitter % This points to the pathloss file that will actually be used
        frequency_band
        antenna_name
        antenna_type
        tx_height
        electrical_downtilt
        mechanical_downtilt
    end

    methods
        
        function print(obj)
            fprintf(' Sector %d: ',obj.id);
            fprintf('%s %ddB %d°\n',obj.antenna.antenna_type,obj.antenna.mean_antenna_gain,obj.azimuth);
            fprintf('  ');
            obj.scheduler.print;
            fprintf('  UEs: ');
            if ~isempty(obj.users)
                current_last_node = obj.users(1);
                while ~isempty(current_last_node.Next)
                    % Do something
                    fprintf('%d ',current_last_node.Data.id);
                    current_last_node = current_last_node.Next;
                end
                % Do something for the last node
                fprintf('%d ',current_last_node.Data.id);
            end
            fprintf('\n');
            fprintf('  '); obj.RB_grid.print;
        end
        
        % Attachs a user to this eNodeB, first checking that the node is
        % not already in the list. It will update the UE's
        % 'attached_eNodeB' variable, effectively binding the UE to this
        % eNodeB. Remember to also add the user to the scheduler, or it
        % will NOT be served!
        function attachUser(obj,user)
            % If the user list is empty
            if isempty(obj.users)
                obj.users = utils.dlnode(user);
                user.attached_eNodeB = obj.parent_eNodeB;
                user.attached_sector = obj.id;
                obj.attached_UEs = obj.attached_UEs + 1;
                % If there are already some users
            else
                % First check if this user is not already in the list
                current_last_node = obj.users;
                user_already_in = false;
                % While the last node has no Next
                while ~isempty(current_last_node.Next)
                    if current_last_node.Data.id==user.id
                        user_already_in = true;
                    end
                    current_last_node = current_last_node.Next;
                end
                % Process the last node
                if current_last_node.Data.id==user.id
                    user_already_in = true;
                end
                % Now current_last_node is the last node
                % Add the new user after it, if not already in the list
                if ~user_already_in
                    new_node = utils.dlnode(user);
                    new_node.insertAfter(current_last_node);
                    obj.attached_UEs = obj.attached_UEs + 1;
                    user.attached_eNodeB = obj.parent_eNodeB;
                    user.attached_sector = obj.id;
                end
            end
        end
        
        % Deattaches a user from this eNodeB. This function does change
        % the user's 'attached_eNodeB' variable. Remember to delete the
        % user from the scheduler also, or nonexistent users will be
        % scheduled!
        function deattachUser(obj,user)
            % If the user list is empty, do nothing
            if ~isempty(obj.users)
                % Process the node list
                current_last_node = obj.users;
                while ~isempty(current_last_node.Next)
                    % Do something
                    if current_last_node.Data.id==user.id
                        current_last_node.Data.attached_eNodeB = [];
                        % In case we are deleting the head
                        if obj.users.Data.id==user.id
                            obj.users = current_last_node.Next;
                        end
                        current_last_node.disconnect;
                        obj.attached_UEs = obj.attached_UEs - 1;
                        return
                    end
                    current_last_node = current_last_node.Next;
                end
                % Do something for the last node
                if current_last_node.Data.id==user.id
                    current_last_node.Data.attached_eNodeB = [];
                    % In case we are deleting the head
                    if obj.users.Data.id==user.id
                        obj.users = current_last_node.Next;
                    end
                    current_last_node.disconnect;
                    obj.attached_UEs = obj.attached_UEs - 1;
                    return
                end
            end
            
            % Also delete the UE from the scheduler
            obj.scheduler.remove_UE(user.id)
        end
        
        % Queries whether a user is attached
        function is_attached = userIsAttached(obj,user)
            % If the user list is empty, return false
            if ~isempty(obj.users)
                % Process the node list
                current_last_node = obj.users;
                while ~isempty(current_last_node.Next)
                    % Do something
                    if current_last_node.Data.id==user.id
                        is_attached = true;
                        return
                    end
                    current_last_node = current_last_node.Next;
                end
                % Do something for the last node
                if current_last_node.Data.id==user.id
                    is_attached = true;
                    return
                end
                is_attached = false;
            else
                is_attached = false;
            end
        end
        
        % Receives and stores the received feedbacks from the UEs
        function receive_UE_feedback(obj)
            current_node = obj.users;
            max_streams  = 2;
            obj.last_received_feedback.UE_id             = zeros(obj.attached_UEs,1);
            obj.last_received_feedback.tx_mode           = zeros(obj.attached_UEs,1); % From what mode the CQI and RI was taken
            obj.last_received_feedback.nCodewords        = zeros(obj.attached_UEs,1); % Relates to the ACK/NACK report. For the 0-delay case, this and tx_mode could differ (ACK is always delayed)
            obj.last_received_feedback.CQI               = zeros(obj.RB_grid.n_RB,max_streams,obj.attached_UEs);
            obj.last_received_feedback.RI                = zeros(obj.attached_UEs,1);
            obj.last_received_feedback.feedback_received = false(obj.attached_UEs,1);
            UE_feedback_idx = 1;

            for i_=1:obj.attached_UEs
                UE_i  = current_node.Data;
                
                % Fill in list of currently attached UEs
                if i_==1
                    obj.UEs_attached_last_TTI = UE_i;
                else
                    obj.UEs_attached_last_TTI(i_) = UE_i;
                end
                
                % Receive the feedback from each user
                UE_id = UE_i.id;
                feedback_u_ = UE_i.uplink_channel.get_feedback;
                % The first TTI, even with 0 delay there is no feedback, as no ACKs are available
                if ~isempty(feedback_u_)
                    
                    % For the zero delay case, substitute the delayed feedback with a zero-delay CQI and (if applicable), RI feedback.
                    % The rest of the feedback data is delayed 1 TTI (no ACK is possible before the reception is done)
                    if obj.zero_delay_feedback
                        % TX mode relates to how the CQI and RI look like.
                        % nCodewords relates to the ACK/NACK size (not overwritten)
                        feedback_u_.tx_mode = UE_i.feedback.tx_mode;
                        feedback_u_.CQI     = UE_i.feedback.CQI;
                        if (feedback_u_.tx_mode==3) || (feedback_u_.tx_mode==4)
                            feedback_u_.RI = UE_i.feedback.RI;
                        end
                    end

                    % Store feedback traces
                    obj.feedback_trace.store(feedback_u_,...
                        UE_id,...
                        obj.parent_eNodeB.id,...
                        obj.id,...
                        obj.parent_eNodeB.clock.current_TTI);
                    
                    % Store accumulated ACK trace. Done separately because it eases
                    % post-processing. It updates the number of correctly
                    % received bits in the trace
                    obj.sector_trace.store_ACK_report(feedback_u_);
                    
                    % Store the last received feedback for all of the attached
                    % users, as it will be needed by the scheduler.
                    % More refined schedulers may need longer "historical" information
                    obj.last_received_feedback.feedback_received(UE_feedback_idx)               = true;
                    obj.last_received_feedback.UE_id(UE_feedback_idx)                           = UE_id;
                    obj.last_received_feedback.tx_mode(UE_feedback_idx)                         = feedback_u_.tx_mode;
                    obj.last_received_feedback.nCodewords(UE_feedback_idx)                      = feedback_u_.nCodewords;
                    obj.last_received_feedback.CQI(:,1:size(feedback_u_.CQI,1),UE_feedback_idx) = feedback_u_.CQI';
                    if (feedback_u_.tx_mode==3) || (feedback_u_.tx_mode==4)
                        obj.last_received_feedback.RI(UE_feedback_idx) = feedback_u_.RI;
                    end
                else
                    obj.last_received_feedback.feedback_received(UE_feedback_idx) = false;
                end
                
                UE_feedback_idx = UE_feedback_idx + 1;
                current_node = current_node.Next;
            end            
        end
        
        % Schedule users in the RB grid for this sector. Modifies the sent
        % resourceBlockGrid object with the user allocation.
        function schedule_users(obj)
            
            % Reset the power allocation to 0. In this way, if no UEs are attached to the scheduler, no power will be transmitted
            if ~obj.always_on
                obj.RB_grid.power_allocation(:) = 0;
            end
            
            % Continue with scheduling
            obj.scheduler.schedule_users(obj.RB_grid,obj.UEs_attached_last_TTI,obj.last_received_feedback);
            % Store traces
            obj.sector_trace.store_after_scheduling(obj.RB_grid);
        end
    end
end
