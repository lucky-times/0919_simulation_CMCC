classdef eNodebSignaling <handle
% Represents the signaling from the eNodeB to each UE.
%
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       TB_CQI           % CQI used for the transmission of each codeword
       TB_size          % size of the current TB, in bits
       N_used_bits = 0;
       packet_parts = [];
       num_assigned_RBs % Number of assigned RBs 
       assigned_RB_map  % Indicates the RBs of this specific UE
       tx_mode          % Transmission mode used (SISO, tx diversity, spatial multiplexing)
       nLayers          % Number of layers for this transmission
       nCodewords       % How many codewords are being sent
       rv_idx           % Redundancy version index (HARQ) for each codeword
       genie_TB_SINR    % Estimated TB SINR as calculated by the eNodeB
       
       % Other signaling, such as X-layer, could be placed here
   end

   methods
   end
end 
