%test.m
%该文件用于修改可变参数并运行平台产生统计数据。

%将之前所有的图形以及变量删除掉，便于观察仿真结果
close all force;
clc;
clear all
clear global;
clear classes;

SYS_config = NR_load_params;%导入默认仿真参数

%% 在默认参数的基础上对一些可变的兴趣参数进行修改调整
% 一些常用的可变参数如下：
SYS_config.shadow_fading_type         = 'claussen';
SYS_config.compact_results_file       = true;
SYS_config.simulation_time_tti        = 1; %仿真经历的时间段
SYS_config.keep_UEs_still             = true; %UE是否保持静止
SYS_config.UE_speed = 5/3.6;% UE运动速度5km/h,(5/3.6)m/s
SYS_config.TTI_length = 1;% Length of a TTI (subframe), in seconds.

% 波束对准的精确度
SYS_config.beam_accuracy_phi_ue = 0; %UE波束赋型水平角（度）
SYS_config.beam_accuracy_theta_ue = 0; %UE波束赋型垂直角（度）
SYS_config.beam_accuracy_phi_bs = 0; %BS波束赋型水平角（度）
SYS_config.beam_accuracy_theta_bs = 0; %BS波束赋型垂直角（度）
SYS_config.beam_accurancy_type = 'Constant'; %Constant表示对准偏差为精度值常量，Unique表示对准偏差为0-精度值的均匀分布
SYS_config.asynchronization = false; %干扰是否不同步

SYS_config.default_shown_GUI_cells = [];
SYS_config.frequency = 4.9e9; %victim系统的频率

SYS_config.photo = true;
SYS_config.PCe_param = 0; %功控的精确度

SYS_config.isDouble = false;% 是否双系统
SYS_config.interference_type = 0;% 0同频干扰，1邻频干扰，2都有，只在双系统下有效
SYS_config.antenna_mode=0;  % 0表示BS与UE均无波束赋型；1表示BS波束赋型，UE无波束赋型；2表示BS与UE均波束赋型
SYS_config.attatch_mode=3;  % 0表示仅考虑路损；1表示路损加BS端单元增益；2表示路损加UE端单元增益；3表示路损加BS和UE端单元增益。用来确定UE归属小区（撒点）
SYS_config.use_cache = false;% 波束赋形是否使用已经存储好的数据
SYS_config.isSave = false; %是否保存波束赋型产生的数据
SYS_config.antenna_element = [8,16];% UMi UMa原来是8,16
SYS_config.antenna_element_InH = [2,2];% InH 30G原来是4,8
SYS_config.antenna_theta_3db = 6;%度

SYS_config.bandwidth = 100e6; %仿真系统的带宽
SYS_config.UE_tx_power = 10^2.3/1000; %UE的发射功率（单位W，此处为23dBm）
SYS_config.Gama = 1; %功控中的gamma参数
SYS_config.n_UE_served_per_BS = 1; % 每个BS同一时间服务的UE数

SYS_config.antenna.antenna_gain_pattern = 'NRAntennaBeamforming'; %波束赋型的类型（可以在此处增加扩展其他的波束赋型图案）
SYS_config.UE_height = 300;%UAV飞行高度 100 200 300 500
SYS_config.UE_r = 1200; %UAV距中心基站的极径  100 300 500 800 1200
str = sprintf("UAV 飞行高度h=%d, 极径r=%d", SYS_config.UE_height, SYS_config.UE_r);
disp(str)
% SYS_config.eNodeB_pos = zeros(19,2);
SYS_config.UE_max_antenna_gain = 5; %UE的最大天线增益（单元增益）
SYS_config.AntennaPattern3d = true;
SYS_config.UE_receiver_noise_figure = 10; %UE端接收机的噪声系数
SYS_config.BS_receiver_noise_figure = 10; %BS端接收机的噪声系数
SYS_config.macroscopic_pathloss_model = 'TS36942'; %第一个系统路损模型的选择，TS38900表示NR系统，TS36942表示LTE路损
SYS_config.macroscopic_pathloss_model2 = 'TS36942'; %第二个系统路损模型的选择，TS38900表示NR系统，TS36942表示LTE路损

SYS_config.asynchronization_switch=0;%判断干扰类型，0表示同步，1表示异步，2表示各50%

SYS_config.cable_loss=0; %天线馈线损耗
SYS_config.scene_type = 'UMA'; %场景参数的选择
SYS_config.beam_loss=0; %波束损耗
SYS_config.tilt = 0; %基站天线上倾角
SYS_config.labelMode = 'Chinese'; %Chinese English 画图模式，图例、坐标轴、标题是采用英文还是中文

switch SYS_config.scene_type
    case 'UMA' %UMA场景下一些参数的设定
        SYS_config.isWraparound = false; %是否使用Wraparound
        SYS_config.isNewUMA = true;% 关系到室内外用户比例
        if SYS_config.isNewUMA % 选定ISD，新版本200 ISD的室内外用户比例为0：10；老版本500ISD的室内外用户比例为8:2
            SYS_config.ISD = 1000; 
        else
            SYS_config.ISD = 500;
        end
        SYS_config.shift_mode = 0;%是否共址,0代表共址，1代表偏移100%,2自己设置偏移量
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %偏置角度，单位度
            SYS_config.ISD_between_systems = (2^0.5)/2*SYS_config.ISD; %偏置间隔，单位米
        end
        SYS_config.hand_over = false; %是否采用3dBhandover
        SYS_config.sector_azimuths = 0:360/3:359; %扇区天线的方向角
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro'; %第一个系统的环境
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro'; %第二个系统的环境
        if SYS_config.ISD <= 300
            if  SYS_config.isWraparound
                SYS_config.map_resolution = 10; %设置宏路损图像素点的粒度
                SYS_config.shadow_fading_map_resolution = 10;%设置阴影衰落图谱像素点的粒度
            else
                SYS_config.map_resolution = 5;
                SYS_config.shadow_fading_map_resolution=5;
            end
        else
            if  SYS_config.isWraparound
                SYS_config.map_resolution = 20;
                SYS_config.shadow_fading_map_resolution = 20;
            else
                SYS_config.map_resolution = 50;
                SYS_config.shadow_fading_map_resolution=10;
            end
        end
        SYS_config.n_snapshot = 1;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 4; %LOS下阴影衰落标准差
        SYS_config.shadow_fading_sd_NLOS = 6; %NLOS下阴影衰落标准差
        SYS_config.r_eNodeBs = .5; %阴影衰落站间相关系数
        SYS_config.site_height=25; %基站的站址高度
        SYS_config.eNodeB_tx_power = 10^4.3/1000; %基站的最大发射功率 w
        SYS_config.antenna.max_antenna_gain = 8; %BS端最大天线增益
        SYS_config.default_shown_GUI_cells = 1:57;
    case 'UMI' %包括Random drop 和manhattan两种子场景
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        else
            SYS_config.map_resolution = 10;
            SYS_config.shadow_fading_map_resolution=10;
        end
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.shift_mode = 0;%是否共址,1代表不共址
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
        SYS_config.n_snapshot = 5;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_NLOS = 7.82;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=10;
        SYS_config.eNodeB_tx_power = 10^3/1000;% 单位w
        SYS_config.isManhattan = false; %是否采用曼哈顿拓扑，true的话采用曼哈顿拓扑，否则采用Random drop场景
        if SYS_config.isManhattan
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = [];
    case 'InH' %包括Open office场景与Home femto场景
        SYS_config.shift_mode = 0;%是否共址,0代表共址
        SYS_config.isFemto = false; %是否为femto场景，若为true则为Home femto场景，否则为Open office场景
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 10;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE的高度为1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% 单位w     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
    case 'RMa'
        SYS_config.isWraparound = true; %是否使用Wraparound
        SYS_config.frequency = 7e9;
        SYS_config.frequency2 = 7e9;
        SYS_config.ISD = 1732;
        SYS_config.shift_mode = 0;%是否共址,0代表共址，1代表偏移100%,2自己设置偏移量
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45;
            SYS_config.ISD_between_systems = 50;
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0:360/3:359;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'rural_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'rural_macro';
        if SYS_config.isWraparound
            SYS_config.map_resolution = 80;
            SYS_config.shadow_fading_map_resolution = 80;
        else
            SYS_config.map_resolution = 50;
            SYS_config.shadow_fading_map_resolution = 50;
        end
        SYS_config.n_snapshot = 10;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_LOS2 = 6; %900中RMA场景下LOS状态有两个LOS。
        SYS_config.shadow_fading_sd_NLOS = 8;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=35;
        SYS_config.eNodeB_tx_power = 10^4.3/1000;% 单位w
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = 1:57;
    case 'UMa_to_UMi'
        SYS_config.isDouble = true;% 是否双系统
        SYS_config.isNewUMA = true;
        SYS_config.isManhattan = false; %微小区是否是曼哈顿，支持UMa与UMi两种子场景的共存研究
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        else
            SYS_config.map_resolution = 10;
            SYS_config.shadow_fading_map_resolution=10;
        end
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.shift_mode = 1;%是否共址,1代表不共址
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% 全向天线
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
        SYS_config.n_snapshot = 5;%设置快照数
        SYS_config.n_snapshot2 = 5;%设置第二个系统快照数
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_NLOS = 7.82;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=10;%10
        SYS_config.eNodeB_tx_power = 10^3.3/1000;% 单位w
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = [];
        
        if SYS_config.isManhattan
            SYS_config.shift_mode = 1;% 两个系统间的偏移，1代表默认（曼哈顿位于中间），2允许设置（第二个系统相对于曼哈顿中心的位置）
            if SYS_config.shift_mode == 2
                SYS_config.angle_between_systems = 45; %偏置角度
                SYS_config.ISD_between_systems = 200; %偏置间隔
            end
        end
        
        SYS_config.sector_azimuths2 = 0:360/3:359; %异构干扰系统的扇区方向角
        SYS_config.eNodeB_tx_power2 = 10^4.3/1000; %异构干扰系统基站的发射功率
        SYS_config.site_height2 = 25;              %异构干扰系统基站的高度
        SYS_config.shadow_fading_sd_LOS_het = 4;   %异构干扰系统基站LOS下的阴影衰落标准差
        SYS_config.shadow_fading_sd_NLOS_het = 6;  %异构干扰系统基站NLOS下的阴影衰落标准差
        SYS_config.antenna.max_antenna_gain2 = 8;  %异构干扰系统基站天线最大单元增益
    case 'UMa_to_InH'
        SYS_config.isDouble = true;% 是否双系统
        SYS_config.isFemto = true;% 支持UMa与InH两种子场景的异构共存研究
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 2;
            SYS_config.shadow_fading_map_resolution=2;
        else
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.shift_mode = 1;% 两个系统间的偏移，1代表默认（inh位于中间），2允许设置（第二个系统相对于inh中心的位置）
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %偏置角度
            SYS_config.ISD_between_systems = 30; %偏置间隔
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% 全向天线
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
        SYS_config.n_snapshot = 5;%设置快照数
        SYS_config.n_snapshot2 = 5;%设置第二个系统快照数
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=3;
        SYS_config.UE_height=1.5;
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% 单位w
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
        SYS_config.sector_azimuths2 = 0:360/3:359;
        SYS_config.eNodeB_tx_power2 = 10^4.3/1000;
        SYS_config.site_height2 = 25;
        SYS_config.shadow_fading_sd_LOS_het = 4;
        SYS_config.shadow_fading_sd_NLOS_het = 6;
        SYS_config.antenna.max_antenna_gain2 = 8;
    case 'UMi_to_InH'
        SYS_config.isDouble = true;% 是否双系统
        SYS_config.isFemto = false;% 仅支持Random drop对InH两种子场景的异构共存研究
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 2;
            SYS_config.shadow_fading_map_resolution=2;
        else
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.shift_mode = 1;% 两个系统间的偏移，1代表默认（inh位于中间），2允许设置（第二个系统相对于inh中心的位置）
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %偏置角度
            SYS_config.ISD_between_systems = 30; %偏置间隔
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% 全向天线
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
        SYS_config.n_snapshot = 5;%设置快照数
        SYS_config.n_snapshot2 = 5;%设置第二个系统快照数
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=3;
        SYS_config.UE_height=1.5;
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% 单位w
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.sector_azimuths2 = 0;
        SYS_config.eNodeB_tx_power2 = 10^3.3/1000;
        SYS_config.site_height2 = 10;
        SYS_config.shadow_fading_sd_LOS_het = 4;
        SYS_config.shadow_fading_sd_NLOS_het = 7.82;
        SYS_config.antenna.max_antenna_gain2 = 8;
        
    case 'InH2' %包括Open office场景与Home femto场景
        SYS_config.shift_mode = 0;%是否共址,0代表共址
        SYS_config.isFemto = false; %是否为femto场景，若为true则为Home femto场景，否则为Open office场景
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 10;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE的高度为1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% 单位w     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
    case 'InH3' %包括Open office场景与Home femto场景
        SYS_config.InH3_d = 10;
        SYS_config.shift_mode = 0;%是否共址,0代表共址
        SYS_config.isFemto = false; %是否为femto场景，若为true则为Home femto场景，否则为Open office场景
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 50;%设置快照数
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;macroscopic_pathloss_model_settings
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE的高度为1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% 单位w     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
end

output_results_file = NR_sim_main(SYS_config); %运行主程序
%simulation_data                   = load(output_results_file); %导入数据
close all;

% GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data); %绘制拓扑图