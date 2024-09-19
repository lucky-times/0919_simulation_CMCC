classdef sectorTrace < handle
% This class stores, for each eNodeB's sector the traces that we wanto to store. eg, CQI assignments,
% throughput, etc, etc.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % whose slot is each
       user_allocation
       % how much power to allocate to each slot, in Watts
       power_allocation
       % Current insertion index. One trace is stored every TTI, so no need
       % to care about "holes". This is equal to the TTI
       TTI_idx
       % Sent data this TTI (bits)
       sent_data
       % Data that was acknowledged (bits)
       acknowledged_data
       % Equivalent but for TBs
       expected_ACKs % Total number of TBs sent (eg. sum of ACKs and NACKs)
       received_ACKs % The ones that were ACKs
   end

   methods
       function obj = sectorTrace(RB_grid_object,maxStreams,simulation_length_TTI)
           %time_slots_power_allocation = size(RB_grid_object.power_allocation,1);
           n_RB = RB_grid_object.n_RB;
           obj.user_allocation   = zeros(n_RB,simulation_length_TTI,'uint16');
           %obj.power_allocation  = zeros(time_slots_power_allocation,n_RB,maxStreams,simulation_length_TTI);
           obj.sent_data         = zeros(maxStreams,simulation_length_TTI,'uint32');
           obj.acknowledged_data = zeros(maxStreams,simulation_length_TTI,'uint32');
           obj.expected_ACKs     = zeros(maxStreams,simulation_length_TTI,'uint16');
           obj.received_ACKs     = zeros(maxStreams,simulation_length_TTI,'uint16');
           obj.TTI_idx = 1;
       end
       % Stores traces after the scheduling was done (Rb grid size)
       function store_after_scheduling(obj,the_RB_grid)
           obj.user_allocation(:,obj.TTI_idx) = the_RB_grid.user_allocation;
           nCodewords = length(the_RB_grid.size_bits);
           obj.sent_data(1:nCodewords,obj.TTI_idx) = the_RB_grid.size_bits;
           obj.TTI_idx = 1 + obj.TTI_idx;
       end
       % Add the trace from a feedback report
       function store_ACK_report(obj,feedback)
           nCodewords   = feedback.nCodewords;
           UE_scheduled = feedback.UE_scheduled;
           ACK          = feedback.ACK;
           TB_size      = feedback.TB_size;
           trace_TTI    = feedback.TTI_idx;
           % In case the ACK report is incorrectly filled. The most important check is whether the user was scheduled or not.
           checked_TB_size = uint32(UE_scheduled .* ACK(:) .* TB_size(:));
           obj.acknowledged_data(:,trace_TTI) = obj.acknowledged_data(:,trace_TTI) + checked_TB_size;
           max_streams = size(obj.expected_ACKs,1);
           if UE_scheduled
               filler_ACK_zeros = zeros(max_streams-nCodewords,1,'uint16');
               current_expected_ACKs = [ones(length(ACK),1,'uint16'); filler_ACK_zeros];
               current_ACKs          = [uint16(ACK(:)); filler_ACK_zeros];
               obj.expected_ACKs(:,trace_TTI) = obj.expected_ACKs(:,trace_TTI) + current_expected_ACKs(:);
               obj.received_ACKs(:,trace_TTI) = obj.received_ACKs(:,trace_TTI) + current_ACKs(:);
           end
       end
   end
end 
