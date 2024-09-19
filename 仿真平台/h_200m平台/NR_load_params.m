function SYS_config = NR_load_params
%% Debug options
SYS_config.debug_level = 1;
SYS_config.compact_results_file = true;
SYS_config.default_shown_GUI_cells = [];

%% Plotting options
SYS_config.photo = true;

%% ACIR设置
SYS_config.ACIR_dB = 15;
SYS_config.isACIR_loop = true;
SYS_config.loop_step = 5;

%% 场景变量
SYS_config.scene_type = 'UMA';
SYS_config.isDouble = true;% 是否双系统
SYS_config.interference_type = 2;% 0同频干扰，1邻频干扰，2都有，只在双系统下有效
SYS_config.antenna_mode=1;  % 0表示BS与UE均无波束赋型；1表示BS波束赋型，UE无波束赋型；2表示BS与UE均波束赋型
SYS_config.attatch_mode=3;  % 1表示路损加BS端单元增益；2表示路损加UE端单元增益；3表示路损加BS和UE端单元增益。用来确定UE归属小区（撒点）
SYS_config.macroscopic_pathloss_model = 'TS38900';
SYS_config.macroscopic_pathloss_is_model = true;
SYS_config.shadow_fading_sd_NLOS_LTE = 10;
SYS_config.UE_per_eNodeB = 5;

SYS_config.isNewUMA = false;
SYS_config.isManhattan = false;
SYS_config.isFemto = false;
SYS_config.isWraparound = false;

SYS_config.nr_eNodeB_rings = 2;

SYS_config.beam_loss=0;
SYS_config.cable_loss=0;
%% 波束赋形设置
SYS_config.use_cache = true;% 波束赋形是否使用已经存储好的数据
SYS_config.isSave = true;
SYS_config.antenna_element = [8,16];% UMi UMa原来是8,16
SYS_config.antenna_element_InH = [4,8];% InH 30G原来是4,8

%% General options
SYS_config.frequency = 26e9;
SYS_config.frequency2 = 262e8;
SYS_config.bandwidth = 200e6;
SYS_config.always_on = true; % Controls whether the eNodeB is always radiating power (dafault and worse-case scenario) or no power is used when no UEs are attached

%% 控制全局随机数的变量
SYS_config.seedRandStream = true;
SYS_config.RandStreamSeed = 0;

%% Simulation time
SYS_config.simulation_time_tti = 1;

%% 阴影衰落
SYS_config.shadow_fading_type = 'claussen'; % Right now only 2D space-correlated shadow fading maps implemented

% Configure the network source
switch SYS_config.shadow_fading_type
    case 'claussen'
        SYS_config.shadow_fading_map_resolution = 5; % Recommended value for 8 neighbors
        SYS_config.shadow_fading_n_neighbors    = 8; % Either 4 or 8
        SYS_config.shadow_fading_mean           = 0;
        SYS_config.shadow_fading_sd             = 10;
        SYS_config.r_eNodeBs                    = 0.5; % inter-site shadow fading correlation
    case 'none'
        % do not use shadow fading
    otherwise
        error([SYS_config.shadow_fading_type ' shadow fading type not supported']);
end

%% 噪声设置
SYS_config.UE.receiver_noise_figure = 11;
SYS_config.UE.BS_receiver_noise_figure = 11;
SYS_config.UE.thermal_noise_density = -174;

%% UE设置
SYS_config.UE_speed = 5/3.6;% 5km/h,(5/3.6)m/s
SYS_config.keep_UEs_still = true;
SYS_config.PCe_param = 0;
SYS_config.Gama = 1;
SYS_config.UE_tx_power = 10^2.24/1000;
SYS_config.UE_max_antenna_gain = 3;

%% eNodeB设置
SYS_config.antenna.antenna_gain_pattern = 'omnidirectionalAntennaBeamforming';
SYS_config.AntennaPattern3d = true;

%% 波束对准的精确度
SYS_config.beam_accuracy_phi_ue=10; %UE波束赋型水平角（度）
SYS_config.beam_accuracy_theta_ue=5; %UE波束赋型垂直角（度）
SYS_config.beam_accuracy_phi_bs=10; %BS波束赋型水平角（度）
SYS_config.beam_accuracy_theta_bs=5; %BS波束赋型垂直角（度）
SYS_config.beam_accurancy_type='Constant'; %Constant表示对准偏差为精度值常量，Unique表示对准偏差为0-精度值的均匀分布

%% Where to save the results
SYS_config.results_folder = './results';
SYS_config.beam_cache_folder = '.\beam_cache\';
SYS_config.results_file = 'auto';

%% Values that should not be changed
% SYS_config.antenna_azimuth_offsett = 30;
SYS_config.antenna_azimuth_offsett = 0;
%% 下面是以前的LTE_load_params_dependant(SYS_config)函数
matlab_release                 = version('-release');
matlab_year                    = str2double(matlab_release(1:(end-1)));
SYS_config.matlab_release_year = matlab_year;

%% Random number generation
if SYS_config.seedRandStream
    if SYS_config.matlab_release_year >= 2011
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',SYS_config.RandStreamSeed));
    else
        RandStream.setDefaultStream(RandStream('mt19937ar','Seed',SYS_config.RandStreamSeed));
    end
else
    if SYS_config.matlab_release_year >= 2011
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',rand*intmax('uint32')));
    else
        RandStream.setDefaultStream(RandStream('mt19937ar','Seed',rand*intmax('uint32')));
    end
end

%% 基站是否总是开启
SYS_config.always_on = true;

%% Ratio of the whole power dedicated to signaling
SYS_config.signaling_ratio = 0; % 0~1

%% Moved from the "do not touch" section
SYS_config.RB_bandwidth = 200e3;         % Frequency in Hz
SYS_config.TTI_length = 5;          % Length of a TTI (subframe), in seconds.

%% Transmission parameters (used for the throughput calculation)
switch SYS_config.bandwidth
    case 1.4e6
        SYS_config.N_RB = 6;
        SYS_config.fft_points = 128;
    case 3e6
        SYS_config.N_RB = 15;
        SYS_config.fft_points = 256;
    case 5e6
        SYS_config.N_RB = 25;
        SYS_config.fft_points = 512;
    case 10e6
        SYS_config.N_RB = 50;
        SYS_config.fft_points = 1024;
    case 15e6
        SYS_config.N_RB = 75;
        SYS_config.fft_points = 1536;
    case 20e6
        SYS_config.N_RB = 100;
        SYS_config.fft_points = 2048;
    case 200e6
        SYS_config.N_RB = 1000;
        SYS_config.fft_points = 20480;
    case 400e6
        SYS_config.N_RB = 1000*2;
        SYS_config.fft_points = 20480*2;
    otherwise
        error('Bandwidth not supported');
end

end