function [eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkShadowFadingMap] = NR_generate_network(SYS_config)
% 网络拓扑产生文件
% 输入参数：
% SYS_config：参数配置
% 输出参数：
% eNodeB_sites：eNodeB实体集
% eNodeB_sectors：扇区实体集
% networkPathlossMap：路损映射矩阵
% networkShadowFadingMap：阴影衰落映射矩阵
%

if SYS_config.debug_level>=1 %输出调试信息
    sprintf('Generating network\n');
end
%% 建立拓扑合配置场景
[eNodeB_sites,eNodeB_sectors,networkPathlossMap] = NR_init_network(SYS_config);

%% 配置导频信号与数据信号的比例，本平台不考虑导频信号
for b_=1:length(eNodeB_sites)
    for s_=1:length(eNodeB_sites(b_).sectors)
        data_power      = eNodeB_sites(b_).sectors(s_).max_power * (1-SYS_config.signaling_ratio);%signaling_ration=0
        signaling_power = eNodeB_sites(b_).sectors(s_).max_power * SYS_config.signaling_ratio;
        eNodeB_sites(b_).sectors(s_).max_power       = data_power;% 发送数据功率，这里就只用这个
        eNodeB_sites(b_).sectors(s_).signaling_power = signaling_power;% 导频功率，这里为0
    end
end

%% 求同频和邻频基站
num_hetnet_sites = networkPathlossMap.num_hetnet_sites;% 异构基站数
num_first_sites = networkPathlossMap.num_first_sites;% 第一个系统基站数
for b_ = 1:length(eNodeB_sites)
    NR_init_get_eNodeB_neighbors(SYS_config,eNodeB_sites(b_), eNodeB_sites,num_hetnet_sites,num_first_sites);
end

%% 产生阴影衰落
networkShadowFadingMap = NR_generate_shadowfading(networkPathlossMap,eNodeB_sites,SYS_config);

%% 给LTE系统加上MCL
% 5G根据提案没有MCL
if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942')
    % 假如第一个系统为LTE
    if SYS_config.isDouble
        num_first_sectors = networkPathlossMap.num_first_sectors;
    else
        num_first_sectors = length(eNodeB_sectors);
    end
    LTE_sectors = 1:num_first_sectors;
    if strcmp(SYS_config.macroscopic_pathloss_model_settings.environment,'urban_macro')
        MCL = 70;% urban 场景MCL为70dB
    elseif strcmp(SYS_config.macroscopic_pathloss_model_settings.environment,'rural_macro')
        MCL = 80;% rural 场景MCL为80dB
    end
%     networkPathlossMap.apply_MCL(LTE_sectors,MCL);
end
if SYS_config.isDouble
    if strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
        % 假如第二个系统为LTE
        num_first_sectors = networkPathlossMap.num_first_sectors;
        LTE_sectors = num_first_sectors+1:length(eNodeB_sectors);
        if strcmp(SYS_config.macroscopic_pathloss_model_settings2.environment,'urban_macro')
            MCL = 70;
        elseif strcmp(SYS_config.macroscopic_pathloss_model_settings2.environment,'rural_macro')
            MCL = 80;
        end
        networkPathlossMap.apply_MCL(LTE_sectors,MCL);
    end
end

%% 基于链路损耗与发射功率求基站覆盖范围
% NR_calculate_cell_coverage(SYS_config,networkPathlossMap,eNodeB_sites,eNodeB_sectors,networkShadowFadingMap);

% To avoid error of missing return argument
if ~exist('networkShadowFadingMap','var')
    networkShadowFadingMap = [];
end
