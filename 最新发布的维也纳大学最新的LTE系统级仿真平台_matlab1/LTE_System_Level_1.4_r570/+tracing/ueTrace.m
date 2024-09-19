classdef ueTrace < handle
% This class stores, for each UE the traces that we wanto to store. eg, CQI,
% throughput, etc, etc.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       latency_time_scale % This is used to compute the averate throughput using an exponential filter
       av_temp 
       TTI_length_s       % Length of one TTI
       nCodewords         % Number of codewords used
       CQI_sent           % CQI thta was fed back
       attached_eNodeB    % Index of the attached eNodeB
       attached_sector    % Index of the attached sector
       position           % UE position
       assigned_RBs       % Number of assigned RBs
       ACK                % ACK/NACK
       TB_CQI             % TB CQI
       TB_size            % TB size
       N_used_bits        % Number of used bits
       BLER               % The used BLER for this TB
       avg_throughput     % Average throughput calculated with an exponential averaging window
       TB_SINR_dB         % The AWGN-equivalent TB SINR as measured by the link quality model. Stored in dB
       SINR_RS_not        % Similar in concept to the RS SINR. Stored in dB
       
       % Optional traces (can be turned on/off in the config file)
       trace_SINR = false;
       SINR            % Stores the RS SINRs [dB]
       SNR             % Stores the average prequalization SNR
   end

   methods
       function obj = ueTrace(simulation_length_TTI,n_RB,maxStreams,traces_config,latency_time_scale,TTI_length_s)
           
           % Trace of sent CQIs is initialized to -1 (an invalid value).
           % Storing more than 1 CQI feedback per TTI is supported.
           if traces_config.unquantized_CQI_feedback
               obj.CQI_sent    = zeros(maxStreams,n_RB,simulation_length_TTI,'single')-1;
           else
               obj.CQI_sent    = zeros(maxStreams,n_RB,simulation_length_TTI,'int8')-1;
           end
           obj.latency_time_scale = latency_time_scale;
           obj.TTI_length_s       = TTI_length_s;
           obj.nCodewords         = zeros(1,simulation_length_TTI,'uint8');
           obj.attached_eNodeB    = zeros(1,simulation_length_TTI,'uint16');
           obj.attached_sector    = zeros(1,simulation_length_TTI,'uint8');
           obj.position           = NaN(2,simulation_length_TTI);
           obj.assigned_RBs       = zeros(1,simulation_length_TTI,'uint8');
           obj.ACK                = false(maxStreams,simulation_length_TTI);
           obj.TB_CQI             = NaN(maxStreams,simulation_length_TTI,'single');
           obj.TB_size            = zeros(maxStreams,simulation_length_TTI,'uint32'); % Max is 80*6*200 (20 MHz, 1 ms for 1 user)
           obj.N_used_bits        = zeros(maxStreams,simulation_length_TTI,'uint32'); % Max is 80*6*200 (20 MHz, 1 ms for 1 user)
           obj.BLER               = zeros(maxStreams,simulation_length_TTI);
           obj.TB_SINR_dB         = zeros(maxStreams,simulation_length_TTI);
           obj.avg_throughput     = zeros(maxStreams,simulation_length_TTI);
           obj.SINR_RS_not        = zeros(1,simulation_length_TTI);
           
           obj.trace_SINR         = false;%traces_config.trace_SINR;
           
       end
       % Trace this specific TTI
       function store(obj,nCodewords,CQI,attached_eNodeB,attached_sector,position,tti_idx,assigned_RBs,ACK,TB_CQI,TB_size,BLER,TB_SINR_dB,N_used_bits,SINR_RS_not,extra_traces)
           % Optional varargin variables to trace are:
           %  - extra_traces{1} -> SINR
           %  - extra_traces{2} -> SNR
           obj.nCodewords(tti_idx)                          = nCodewords;
           obj.CQI_sent(1:size(CQI,1),:,tti_idx)            = CQI;
           obj.attached_eNodeB(tti_idx)                     = attached_eNodeB;
           obj.attached_sector(tti_idx)                     = attached_sector;
           obj.position(:,tti_idx)                          = position;
           obj.assigned_RBs(tti_idx)                        = assigned_RBs;
           obj.ACK(1:length(ACK),tti_idx)                   = ACK;
           obj.TB_CQI(1:length(TB_CQI),tti_idx)             = TB_CQI;
           obj.TB_size(1:length(TB_size),tti_idx)           = TB_size;
           obj.N_used_bits(1:length(N_used_bits),tti_idx)   = N_used_bits;
           obj.BLER(1:length(BLER),tti_idx)                 = BLER;
           obj.TB_SINR_dB(1:length(TB_SINR_dB),tti_idx)     = TB_SINR_dB;
           obj.SINR_RS_not(tti_idx)                         = SINR_RS_not;
           
           throughput     = zeros(size(obj.avg_throughput,1),1);
           throughput(1:length(ACK),1) = TB_size(:).*ACK(:) / obj.TTI_length_s;

           obj.av_temp = min(obj.latency_time_scale,tti_idx);
           if tti_idx==1
               obj.avg_throughput(1:length(throughput),tti_idx) = throughput;
           else
               obj.avg_throughput(1:length(throughput),tti_idx) = (1-1/obj.av_temp)*obj.avg_throughput(1:length(throughput),tti_idx-1) + 1/obj.av_temp*throughput(:);
           end
           if obj.trace_SINR
               obj.SINR(:,:,tti_idx) = extra_traces{1}; % In dB
               obj.SNR(tti_idx)      = extra_traces{2}; % In dB
           end
       end
   end
end 
