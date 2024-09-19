classdef proportionalFairScheduler < schedulers.lteScheduler
% A proportional fair LTE scheduler.
% (c) Josep Colom Ikuno, INTHFT, 2010

   properties
       % See the lteScheduler class for a list of inherited attributes
       alpha
       beta
   end

   methods
       
       % Class constructor.
       % The needed parameters are the following:
       %   - window_size: in order to calculate the eNodeB's past throughput, a time window is used. Specified in TTIs
       %   - alpha: \alpha parameter that is used to calculate the priority P
       %   - beta:  \beta parameter that is used to calculate the priority P
       function obj = proportionalFairScheduler(scheduler_params,attached_eNodeB_sector,varargin)
           % Fill in basic parameters (handled by the superclass constructor)
           obj      = obj@schedulers.lteScheduler(scheduler_params,attached_eNodeB_sector);
           obj.name = 'Proportional fair scheduler';
           if isempty(varargin)
               obj.alpha = 1;
               obj.beta  = 1;
           elseif length(varargin)==2
               obj.alpha = varargin{1};
               obj.beta  = varargin{2};
           else
               error('Valid constructors are: proportionalFairScheduler(attached_eNodeB_sector)->alpha=1,beta=1 or proportionalFairScheduler(attached_eNodeB_sector,alpha,beta)');
           end
       end
       
       % Dummy functions required by the lteScheduler Abstract class implementation
       % Add UE (no UE memory, so empty)
       function add_UE(obj,UE_id)
       end
       % Delete UE (no UE memory, so empty)
       function remove_UE(obj,UE_id)
       end
       
       % Schedule the users in the given RB grid
       function schedule_users(obj,RB_grid,attached_UEs,last_received_feedbacks)
           % Power allocation
           % Nothing here. Leave the default one (homogeneous)
           
           RB_grid.size_bits = 0;
           
           % For now use the static tx_mode assignment
           RB_grid.size_bits = 0;
           nCodewords  = RB_grid.nCodewords;
           nLayers     = RB_grid.nLayers;
           tx_mode     = RB_grid.tx_mode;
           current_TTI = obj.clock.current_TTI;

           if ~isempty(attached_UEs)
               
               % Fill in P matrix
               P = zeros(size(last_received_feedbacks.CQI,1),size(last_received_feedbacks.CQI,2));
               CQI_efficiency = obj.get_spectral_efficiency(last_received_feedbacks.CQI);
               R_k_RB = sum(180e3*CQI_efficiency/obj.clock.TTI_time,3);  % Requested throughput for all of the spatial streams
               TTI_to_read = max(current_TTI-obj.feedback_delay_TTIs,1); % Realistically read the ACKed throughput
               T_k = zeros(length(attached_UEs),1);
               for u_=1:length(attached_UEs)
                   UE_id = last_received_feedbacks.UE_id(u_);
                   if UE_id~=0
                       T_k(u_) = sum(obj.UE_traces(UE_id).avg_throughput(:,TTI_to_read)); % Mean throughput, averaged with an exponential window)
                   else
                       T_k(u_) = Inf;
                   end
               end
               T_k_mat = repmat(T_k,[1 RB_grid.n_RB]);
               % Substitute 0 average throughput by the minimum allowable
               % value to avoid NaNs when R_k==0
               T_k_mat(T_k_mat==0) = eps;
               
               P = R_k_RB.^obj.alpha ./ T_k_mat.^obj.beta;
               if current_TTI ~= 1 % at the first TTI there is no feedback --> it does not work, but it is necessary for zero delay feedback
                    UE_id_list = obj.get_max_UEs(P,last_received_feedbacks.UE_id);
               else
                    UE_id_list = obj.get_max_UEs(ones(size(P)),1:length(attached_UEs));
               end
               
               % Fill in RB grid
               RB_grid.user_allocation(:) = UE_id_list;

               % CQI assignment. TODO: implement HARQ
               obj.schedule_users_common(RB_grid,attached_UEs,last_received_feedbacks,current_TTI,nLayers,nCodewords,tx_mode);
           end
       end
       
       function schedule_users_common(obj,RB_grid,attached_UEs,last_received_feedbacks,current_TTI,nLayers,nCodewords,tx_mode)
           schedule_users_common@network_elements.lteScheduler(obj,RB_grid,attached_UEs,last_received_feedbacks,current_TTI,nLayers,nCodewords,tx_mode);
       end
   end
end 
