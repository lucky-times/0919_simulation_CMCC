function SYS_config = NR_load_params_again( SYS_config )
% 防止某些参数配置出错，这些参数影响到拓扑的建立和波束赋形文件的名字
% 输入参数：
% SYS_config：参数配置
%
% 输出参数：
% SYS_config：参数配置

% check
try
    if SYS_config.n_snapshot~=SYS_config.n_snapshot2 % 暂不支持异构时两种小区不同UE数目的情况
        SYS_config.n_snapshot2 = SYS_config.n_snapshot;
    end
catch
    % 若SYS_config.n_snapshot2没有定义，说明不是异构场景
    % do nothing
end

% 对不同场景的参数进行限制
switch SYS_config.scene_type
    case 'UMA'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;% 快照数与基站同时服务的UE数相乘得到一个小区的UE数
        SYS_config.isManhattan = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case 'UMI'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.shift_mode = 1;
        SYS_config.isNewUMA = false;
        SYS_config.isWraparound = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case {'InH','InH2','InH3'}
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.ISD = [];
        SYS_config.shift_mode = 0;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
        if SYS_config.isFemto
            SYS_config.shift_mode = 1;
        end
        SYS_config.isS2F = true;
    case 'RMa'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.isManhattan = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case 'UMa_to_UMi'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = true;
        SYS_config.isWraparound = false;
        SYS_config.isFemto = false;
        if SYS_config.isManhattan
            SYS_config.ISD = 500;
            SYS_config.map_resolution = 15;
            SYS_config.shadow_fading_map_resolution=15;
        end
    case 'UMa_to_InH'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
    case 'UMi_to_InH'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
end

% 保存的文件命名
name_generator = utils.naming;
SYS_config.results_file = name_generator.results_file(SYS_config);
SYS_config.beam_cache_file = name_generator.beam_cache_file(SYS_config);

end

