%% Automagically filled params and config params that will probably not be changed

% Hardcode the release number
LTE_config.release = 'r427';

%% Check whether the parallel computing toolbox is installed
v = ver;
if isempty(strfind([v.Name],'Parallel Computing Toolbox'))
    LTE_config.parallel_toolbox_installed = false;
else
    LTE_config.parallel_toolbox_installed = true;
end

%% Random number generation
if LTE_config.seedRandStream
    RandStream.setDefaultStream(RandStream('mt19937ar','Seed',LTE_config.RandStreamSeed));
else
    RandStream.setDefaultStream(RandStream('mt19937ar','Seed',rand*intmax('uint32')));
end

%% Backwards-compatibility safeguards

% Minimum coupling loss
if ~isfield(LTE_config,'minimum_coupling_loss')
    LTE_config.minimum_coupling_loss = 0; % Mimic old behavior
end

% Interference
if ~isfield(LTE_config,'always_on')
    LTE_config.always_on = true; % Mimic old behavior
end

% Manual ROI setting
if ~isfield(LTE_config,'manually_set_ROI')
    LTE_config.manually_set_ROI = false; % Mimic old behavior
else
    if LTE_config.manually_set_ROI && (~isfield(LTE_config,'roi_x') || ~isfield(LTE_config,'roi_y'))
        error('When manually setting the ROI you need to specify the roi_x and roi_y variables');
    end
end

% UE antenna gain (assumed omnidirectional)
if ~isfield(LTE_config.UE,'antenna_gain')
    LTE_config.UE.antenna_gain = 0; % Mimic old behavior
end

% Whether to save results in a very compact way or the "traditional way" (results object)
if ~isfield(LTE_config,'compact_results_file') % Save just the UE throughput
    LTE_config.compact_results_file = false; % Mimic old behavior
end

% Traffic map upscaling factor
if ~isfield(LTE_config,'traffic_map_upscaling')
    LTE_config.traffic_map_upscaling = 1; % Mimic old behavior
end

% Whether to delete the pathloss data at the end (saves some space)
if ~isfield(LTE_config,'delete_pathloss_at_end')
    LTE_config.delete_pathloss_at_end = false; % Mimic old behavior
end

% Whether the MCL is re-applied at the end regardless of what is load from
% cache. If set to 'false', this extra step will not be performed when
% loading from cache or planning tool data (slightly shorter startup time)
if ~isfield(LTE_config,'reapply_MCL_to_cache_or_planning_tool')
    LTE_config.reapply_MCL_to_cache_or_planning_tool = false; % Mimic old behavior
end

% Configurable additional penetration loss (e.g. 7dB extra for "in car" losses)
if ~isfield(LTE_config,'additional_penetration_loss') || isempty(LTE_config.additional_penetration_loss)
    LTE_config.additional_penetration_loss = 0; % Mimic old behavior
end
if ~isnumeric(LTE_config.additional_penetration_loss)
    switch LTE_config.additional_penetration_loss
        case 'deep indoor'
            LTE_config.additional_penetration_loss = 23;
        case 'indoor'
            LTE_config.additional_penetration_loss = 17;
        case 'incar'
            LTE_config.additional_penetration_loss = 7;
        case 'outdoor'
            LTE_config.additional_penetration_loss = 0;
        otherwise
            error('"deep indoor", "indor", "incar", and "outdoor" values or a penetration loss value in dB are valid');
    end
end

% Preserve the previously-loaded pregenerated_ff object to avoid continuous
% reading of large files during long batches of short simulations
if ~isfield(LTE_config,'reuse_pregenerated_ff_trace_from_last_run')
    LTE_config.reuse_pregenerated_ff_trace_from_last_run = false; % Mimic old behavior
end

% Suffix to the output results file
if ~isfield(LTE_config,'output_filename_suffix') || isempty(LTE_config.output_filename_suffix)
    suffix_to_use = []; % Mimic old behavior
else
    suffix_to_use = ['_' LTE_config.output_filename_suffix];
end

% Ratio of the whole power dedicated to signaling
if ~isfield(LTE_config,'signaling_ratio')
    LTE_config.signaling_ratio = 0; % Mimic old behavior
else
    if LTE_config.signaling_ratio>1 || LTE_config.signaling_ratio<0
        error('Signaling ratio must be between 0 and 1.');
    end
end

%% Channel trace generation configuration
switch LTE_config.channel_model.type
    case 'winner+'
        winner_model_path = './Winner Channel Model';
        if LTE_config.recalculate_fast_fading || ~exist([LTE_config.pregenerated_ff_file '.mat'],'file')
            channel_gain_wrappers.winnerChannelFactory.prepare_winner_channel_model_path(winner_model_path);
        end
        LTE_config.trace_params = channel_gain_wrappers.winnerChannelFactory.get_default_config(LTE_config.frequency,LTE_config.nTX,LTE_config.nRX,LTE_config.UE_speed);
    otherwise
        LTE_config.trace_params = channel_gain_wrappers.pdpChannelFactory.get_default_config(LTE_config.frequency,LTE_config.nTX,LTE_config.nRX,LTE_config.UE_speed,LTE_config.channel_model.type);
end

%% Moved from the "do not touch" section
LTE_config.RB_bandwidth    = 180e3;         % Frequency in Hz
LTE_config.TTI_length      = 1e-3;          % Length of a TTI (subframe), in seconds.
LTE_config.cyclic_prefix   = 'normal';      % 'normal' or 'extended' cyclic prefix
LTE_config.maxStreams      = 2;             % Maximum number of codewords per TTI. For LTE, that's 2.
                                            % Took the name from HSDPA. In
                                            % the LTE standard is actually
                                            % referred as 'codewords'

%% Warning about 0-delay feedback
if LTE_config.feedback_channel_delay==0
    warning('The assumed allocated power when calculating the SINR and performing the simulation will always be assumed to be %dW. The assigned power will affect the measured SINR, but in order to measure it, the scheduling has to be done first, which cannot be done without proper feedback. Therefore, power assignment is not used for 0-delay.',LTE_config.eNodeB_tx_power);
    LTE_config.always_on = true;
end

%% Unquantized feedback
LTE_config.traces_config.unquantized_CQI_feedback = LTE_config.unquantized_CQI_feedback;

%% Check SISO mode settings
if LTE_config.tx_mode==1 && LTE_config.nTX~=1 && LTE_config.nRX~=1
    error('For SISO mode, nTX and nRX must be set to 1.');
end

%% Check whether the macroscopic pathloss comes from a model or from data
switch LTE_config.macroscopic_pathloss_model
    case 'Capesso'
        LTE_config.macroscopic_pathloss_is_model = false;
    otherwise
        LTE_config.macroscopic_pathloss_is_model = true;
end

%% Date-time string
the_date = clock;
date_time_string = sprintf('%04d%02d%02d_%02d%02d%02d',...
        the_date(1),...                     % Date: year
        the_date(2),...                     % Date: month
        the_date(3),...                     % Date: day
        the_date(4),...                     % Date: hour
        the_date(5),...                     % Date: minutes
        floor(the_date(6)));                % Date: seconds

%% Results file filename
if strcmp(LTE_config.results_file,'auto')
    
    if LTE_config.frequency/1e9 >= 1
        this_freq = sprintf('%3.2fGHz',LTE_config.frequency/1e9);
    else
        this_freq = sprintf('%3.0fMHz',LTE_config.frequency/1e6);
    end
    
    if LTE_config.bandwidth/1e9 >= 1
        this_bw = sprintf('%3.2fGHz',LTE_config.bandwidth/1e9);
    else
        this_bw = sprintf('%3.0fMHz',LTE_config.bandwidth/1e6);
    end
    
    switch LTE_config.tx_mode
        case 1
            MIMO_mode = 'SISO';
        case 2
            MIMO_mode = sprintf('%dx%dTxD',LTE_config.nTX,LTE_config.nRX);
        case 3
            MIMO_mode = sprintf('%dx%dOLSM',LTE_config.nTX,LTE_config.nRX);
        case 4
            MIMO_mode = sprintf('%dx%dCLSM',LTE_config.nTX,LTE_config.nRX);
    end
    
    LTE_config.results_file = fullfile(LTE_config.results_folder,...
        sprintf('%s_freq_%s_bw_%s_%.1fKmph_%dTTIs_%s_%s_%s_%s%s.mat',...
        this_freq,...                             % Frequency
        this_bw,...                               % System Bandwidth
        LTE_config.channel_model.type,...         % Channel model userd
        LTE_config.UE_speed*3.6,...               % Channel speed (Km/h)
        LTE_config.simulation_time_tti,...        % Simulaton length
        date_time_string,...                      % Date string
        strrep(LTE_config.scheduler,' ','_'),...  % Scheduler type
        MIMO_mode,...                             % MIMO mode
        LTE_config.release,...                    % Release number
        suffix_to_use));                          % Optional suffix
    clear suffix_to_use
        
else
     LTE_config.results_file = fullfile(LTE_config.results_folder,[LTE_config.results_file '.mat']);
end

%% Macroscopic pathloss cache filename
if strcmp(LTE_config.antenna.antenna_gain_pattern, 'kathreinTSAntenna')
    LTE_config.AntennaPattern3d = true;
end
if strcmp(LTE_config.network_cache,'auto')
    if LTE_config.frequency >= 1e9
        this_freq = sprintf('%3.2fGHz',LTE_config.frequency/1e9);
    else
        this_freq = sprintf('%3.0fMHz',LTE_config.frequency/1e6);
    end
    if strcmp(LTE_config.network_source, 'generated')
        if isfield(LTE_config,'AntennaPattern3d') && LTE_config.AntennaPattern3d
            downtilting = sprintf('_%d°-%d°',LTE_config.antenna.electrical_downtilt,LTE_config.antenna.mechanical_downtilt);
        else
            downtilting = [];
        end
        LTE_config.network_cache = fullfile('./data_files',...
            sprintf('network_%d_rings_%dm_res_%s_%s_antenna%s_%s_freq.mat',...
            LTE_config.nr_eNodeB_rings,...
            LTE_config.map_resolution,...
            strrep(strtrim([LTE_config.macroscopic_pathloss_model ' ' LTE_config.macroscopic_pathloss_model_settings.environment]),' ','_'),...
            LTE_config.antenna.antenna_gain_pattern,...
            downtilting,...
            this_freq));
    else
        % When using planning tool Capesso
        LTE_config.network_cache = fullfile('./data_files',...
            sprintf('network_%dm_res_%s_%s_freq.mat',...
            LTE_config.map_resolution,...
            strrep(strtrim([LTE_config.macroscopic_pathloss_model ' ' LTE_config.macroscopic_pathloss_model_settings.environment]),' ','_'),...
            this_freq));
    end
else
    % Do nothing
end

%% UE position cache filename
if strcmp(LTE_config.UE_cache_file,'auto')
    if ~LTE_config.UEs_only_in_target_sector
        LTE_config.UE_cache_file = fullfile('./data_files',...
            sprintf('UE_cache_%drings_%dUEs_sector_%s.mat',LTE_config.nr_eNodeB_rings,LTE_config.UE_per_eNodeB,date_time_string));
    else
        LTE_config.UE_cache_file = fullfile('./data_files',...
            sprintf('UE_cache_%drings_target_sector_only_%dUEs_sector_%s.mat',LTE_config.nr_eNodeB_rings,LTE_config.UE_per_eNodeB,date_time_string));
    end
else
    % Do nothing
end

%% Fast fading filenames
if strcmp(LTE_config.pregenerated_ff_file,'auto')
    LTE_config.pregenerated_ff_file           = fullfile('./data_files',...
        sprintf('%dx%d_%s_%3.1fMHz_%dKmph_%3.1fs_trace_%s.mat',...
        LTE_config.nTX,...
        LTE_config.nRX,...
        LTE_config.channel_model.type,...
        LTE_config.bandwidth/1e6,...
        LTE_config.UE_speed*3.6,...
        LTE_config.channel_model.trace_length,...
        date_time_string));
end

%% Transmission parameters (used for the throughput calculation)
% We will assume subcarrier spacing of 15 kHz
switch LTE_config.cyclic_prefix
    case 'normal'
        LTE_config.N_sym = 7;
    case 'extended'
        LTE_config.N_sym = 6;
    otherwise
        error('CP can only be "normal" or "extended"');
end
switch LTE_config.bandwidth
    case 1.4e6
        LTE_config.N_RB = 6;
        LTE_config.fft_points = 128;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 9;
            case 'extended'
                LTE_config.CP_length_samples = 32;
        end
    case 3e6
        LTE_config.N_RB = 15;
        LTE_config.fft_points = 256;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 18;
            case 'extended'
                LTE_config.CP_length_samples = 64;
        end
    case 5e6
        LTE_config.N_RB = 25;
        LTE_config.fft_points = 512;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 36;
            case 'extended'
                LTE_config.CP_length_samples = 128;
        end
    case 10e6
        LTE_config.N_RB = 50;
        LTE_config.fft_points = 1024;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 72;
            case 'extended'
                LTE_config.CP_length_samples = 256;
        end
    case 15e6
        LTE_config.N_RB = 75;
        LTE_config.fft_points = 1536;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 108;
            case 'extended'
                LTE_config.CP_length_samples = 384;
        end
    case 20e6
        LTE_config.N_RB = 100;
        LTE_config.fft_points = 2048;
        switch LTE_config.cyclic_prefix
            case 'normal'
                LTE_config.CP_length_samples = 144;
            case 'extended'
                LTE_config.CP_length_samples = 512;
        end
    otherwise
        error('Bandwidth not supported');
end
LTE_config.Ntot = LTE_config.N_RB*12;
LTE_config.fs = 15e3*LTE_config.fft_points;
switch LTE_config.nTX   % number of reference symbols
    case 1
        numb = 4;
    case 2
        numb = 8;
    otherwise
        numb = 12;
end
LTE_config.sym_per_RB_nosync = 12*LTE_config.N_sym - numb; % no synchronization signals 
LTE_config.sym_per_RB_sync = 12*LTE_config.N_sym -(numb+12); % take also synchronization signals into account
LTE_config.scheduler_params.overhead_ref = numb;
LTE_config.scheduler_params.overhead_sync = numb+12;
%% BLER curves
LTE_config.BLER_curves.folder = fullfile(pwd,'data_files','AWGN_BLERs');
LTE_config.BLER_curves.filenames = {
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi1.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi2.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi3.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi4.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi5.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi6.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi7.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi8.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi9.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi10.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi11.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi12.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi13.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi14.mat')
    fullfile(LTE_config.BLER_curves.folder,'AWGN_1.4MHz_SISO_cqi15.mat')
};

%% Store CQI parameters
% CQI 1 is index 1, CQI 2 is index 2, etc...
LTE_config.CQI_params(1).CQI = 1;
LTE_config.CQI_params(1).modulation = 'QPSK';
LTE_config.CQI_params(1).modulation_order = 2;
LTE_config.CQI_params(1).coding_rate_x_1024 = 78;
LTE_config.CQI_params(1).efficiency = 0.1523;

LTE_config.CQI_params(2).CQI = 2;
LTE_config.CQI_params(2).modulation = 'QPSK';
LTE_config.CQI_params(2).modulation_order = 2;
LTE_config.CQI_params(2).coding_rate_x_1024 = 120;
LTE_config.CQI_params(2).efficiency = 0.2344;

LTE_config.CQI_params(3).CQI = 3;
LTE_config.CQI_params(3).modulation = 'QPSK';
LTE_config.CQI_params(3).modulation_order = 2;
LTE_config.CQI_params(3).coding_rate_x_1024 = 193;
LTE_config.CQI_params(3).efficiency = 0.3770;

LTE_config.CQI_params(4).CQI = 4;
LTE_config.CQI_params(4).modulation = 'QPSK';
LTE_config.CQI_params(4).modulation_order = 2;
LTE_config.CQI_params(4).coding_rate_x_1024 = 308;
LTE_config.CQI_params(4).efficiency = 0.6016;

LTE_config.CQI_params(5).CQI = 5;
LTE_config.CQI_params(5).modulation = 'QPSK';
LTE_config.CQI_params(5).modulation_order = 2;
LTE_config.CQI_params(5).coding_rate_x_1024 = 449;
LTE_config.CQI_params(5).efficiency = 0.8770;

LTE_config.CQI_params(6).CQI = 6;
LTE_config.CQI_params(6).modulation = 'QPSK';
LTE_config.CQI_params(6).modulation_order = 2;
LTE_config.CQI_params(6).coding_rate_x_1024 = 602;
LTE_config.CQI_params(6).efficiency = 1.1758;

LTE_config.CQI_params(7).CQI = 7;
LTE_config.CQI_params(7).modulation = '16QAM';
LTE_config.CQI_params(7).modulation_order = 4;
LTE_config.CQI_params(7).coding_rate_x_1024 = 378;
LTE_config.CQI_params(7).efficiency = 1.4766;

LTE_config.CQI_params(8).CQI = 8;
LTE_config.CQI_params(8).modulation = '16QAM';
LTE_config.CQI_params(8).modulation_order = 4;
LTE_config.CQI_params(8).coding_rate_x_1024 = 490;
LTE_config.CQI_params(8).efficiency = 1.9141;

LTE_config.CQI_params(9).CQI = 9;
LTE_config.CQI_params(9).modulation = '16QAM';
LTE_config.CQI_params(9).modulation_order = 4;
LTE_config.CQI_params(9).coding_rate_x_1024 = 616;
LTE_config.CQI_params(9).efficiency = 2.4063;

LTE_config.CQI_params(10).CQI = 10;
LTE_config.CQI_params(10).modulation = '64QAM';
LTE_config.CQI_params(10).modulation_order = 6;
LTE_config.CQI_params(10).coding_rate_x_1024 = 466;
LTE_config.CQI_params(10).efficiency = 2.7305;

LTE_config.CQI_params(11).CQI = 11;
LTE_config.CQI_params(11).modulation = '64QAM';
LTE_config.CQI_params(11).modulation_order = 6;
LTE_config.CQI_params(11).coding_rate_x_1024 = 567;
LTE_config.CQI_params(11).efficiency = 3.3223;

LTE_config.CQI_params(12).CQI = 12;
LTE_config.CQI_params(12).modulation = '64QAM';
LTE_config.CQI_params(12).modulation_order = 6;
LTE_config.CQI_params(12).coding_rate_x_1024 = 666;
LTE_config.CQI_params(12).efficiency = 3.9023;

LTE_config.CQI_params(13).CQI = 13;
LTE_config.CQI_params(13).modulation = '64QAM';
LTE_config.CQI_params(13).modulation_order = 6;
LTE_config.CQI_params(13).coding_rate_x_1024 = 772;
LTE_config.CQI_params(13).efficiency = 4.5234;

LTE_config.CQI_params(14).CQI = 14;
LTE_config.CQI_params(14).modulation = '64QAM';
LTE_config.CQI_params(14).modulation_order = 6;
LTE_config.CQI_params(14).coding_rate_x_1024 = 873;
LTE_config.CQI_params(14).efficiency = 5.1152;

LTE_config.CQI_params(15).CQI = 15;
LTE_config.CQI_params(15).modulation = '64QAM';
LTE_config.CQI_params(15).modulation_order = 6;
LTE_config.CQI_params(15).coding_rate_x_1024 = 948;
LTE_config.CQI_params(15).efficiency = 5.5547;

%% Plot options
LTE_config.plots.BLER_curves                    = 1;
LTE_config.plots.CQI_mapping                    = 2;
LTE_config.plots.antenna_gain_pattern           = 3;
LTE_config.plots.macroscopic_pathloss           = 4;

LTE_config.plots.macroscopic_pathloss_sector(1) = 5;
LTE_config.plots.macroscopic_pathloss_sector(2) = 6;
LTE_config.plots.macroscopic_pathloss_sector(3) = 7;

%LTE_config.plots.sector_assignment              = 10;
%LTE_config.plots.sector_assignment_no_shadowing = 11;
LTE_config.plots.shadow_fading_loss             = 12;
%LTE_config.plots.shadow_fading_loss_histogram   = 13;
LTE_config.plots.initial_UE_positions           = 14;
LTE_config.plots.user_positions                 = 15;
LTE_config.plots.sector_SINR                    = 16;
LTE_config.plots.sector_SINR_no_shadowing       = 17;
LTE_config.plots.sector_SINR_cdf                = 18;
LTE_config.plots.capesso_vs_pathloss_models     = 19;
LTE_config.plots.capesso_maps_begin             = 100;

%% Traffic models
if isfield(LTE_config,'traffic_models')
    if LTE_config.traffic_models.usetraffic_model && ~strcmp(LTE_config.scheduler,'constrained')
        warning('Traffic models are just supported with the constrained scheduler - deactivating traffic models');
        LTE_config.traffic_models.usetraffic_model = false;
    end
else
    LTE_config.traffic_models.usetraffic_model = false;
end

% Some after-cleaning
clear the_date;
clear date_time_string;
clear this_freq;
clear v;