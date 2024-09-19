classdef resourceBlockGrid < handle
% Represents the recource block grid that the scheduler must allocate every
% TTI. Frequency length will depend on the frequency bandwidth. Time length
% will always be 1 ms.
% The grids are organized in the following way:
%
% |<----frequency---->
%  ____ ____ ____ ____  ___
% |____|____|____|____|  |  time (2 subframes)
% |____|____|____|____| _|_
%
% Where the frequency width obviously depends on the allocated bandwidth
% and the time-dimension width is always 2. Access the grid in the
% following way (example for the user allocation):
%
% user_allocation(time_index,frequency_index);
%
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % whose slot is each
       user_allocation
       % how much power to allocate to each slot, in Watts
       power_allocation
       % equivalent, but for the signaling
       power_allocation_signaling
       % number of RB (freq domain)
       n_RB
       % number of symbols per Resource Block (RB), which is 12 subcarriers and 0.5 ms
       sym_per_RB_nosync
       sym_per_RB_sync
       % total size of the RB grid in bits
       size_bits

       tx_mode % Left it here since right now the simulator does not allow for mixed-mode simulations
   end

   methods
       
       % Class constructor and initialisation. Initialise CQIs to 0
       function obj = resourceBlockGrid(n_RB,sym_per_RB_nosync,sym_per_RB_sync)
           max_codewords = 2;
           obj.user_allocation  = zeros(n_RB,1,'uint16'); % We will assume that the streams cannot be scheduled to different UEs.
           obj.power_allocation = zeros(n_RB,1); % TTI-based power allocation. No slot-based power allocation
           obj.power_allocation_signaling = zeros(n_RB,1); % Simplification: equal for each RB
           obj.n_RB              = n_RB;
           obj.sym_per_RB_nosync = sym_per_RB_nosync;
           obj.sym_per_RB_sync   = sym_per_RB_sync;
           obj.size_bits         = zeros(1,max_codewords);
       end
       
       % Sets the power allocation to a default value. Useful for setting a homogeneous power allocation
       function set_homogeneous_power_allocation(obj,power_in_watts_data,power_in_watts_signaling)
           power_per_RB_data                 = power_in_watts_data / obj.n_RB;
           power_per_RB_signaling            = power_in_watts_signaling / obj.n_RB;
           obj.power_allocation(:)           = power_per_RB_data;
           obj.power_allocation_signaling(:) = power_per_RB_signaling;
       end
       function print(obj)
           fprintf('n_RB=%d\n',obj.n_RB);
       end
   end
end 
