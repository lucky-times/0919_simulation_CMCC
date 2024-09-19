function [ UEs ] = NR_generate_users(SYS_config,eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkClock)
% 撒点并创建UE
% 输入参数：
% SYS_config：参数配置
% eNodeB_sites：基站
% eNodeB_sectors：扇区
% networkPathlossMap：路损映射矩阵
% networkClock时钟实体
%
% 输出参数：
% UEs：用户
%

%% 撒点
% UE_spatial_distribution = spatial_distributions.constantUesPerCellUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors,SYS_config);% 实例化UE撒点方式
% UE_positions = UE_spatial_distribution.generate_UE_positions(SYS_config);% UE_positions是所有UE的位置矩阵，这是还没有UE实体

% r = SYS_config.UE_r;
% step = 0.5;
% theta = 0:step:360-step;
UE_positions = [SYS_config.ISD/2; 0]';

%% 创建用户
UEs = network_elements.UE;
UE_pos_vector = zeros(size(UE_positions,1),2);
for u_ = 1:size(UE_positions,1)
    % 用户基本配置
    UEs(u_)     = network_elements.UE; %实体建立
    UEs(u_).id  = u_; %用户id设置
%     UEs(u_).pos = NR_common_pixel_to_pos( UE_positions(u_,:), networkPathlossMap.coordinate_origin, networkPathlossMap.data_res); %用户位置配置
    UEs(u_).pos = UE_positions(u_,:);
    UE_pos_vector(u_,:) = UEs(u_).pos;
    UEs(u_).walking_model = walking_models.straightWalkingModel(SYS_config.UE_speed*SYS_config.TTI_length); % UE运动模型  %与函数模板不匹配，模板参数为方向，此时输入为路程
end
networkPathlossMap.UE_pos_vector = UE_pos_vector;

%% 初始化RB分配
% RB这块从原来的平台继承过来，只用到了使用RB来确定基站给用户分配的带宽
% 如果要加入资源分配，可以适当修改这块
NR_init_RB_grid(SYS_config,eNodeB_sites);

%% UE其他部分的初始化
for u_=1:length(UEs)
    
    % UE天线
    antennas.antenna.attach_antenna_to_UE(UEs(u_),SYS_config);
    
    % 传入噪声系数属性
    UEs(u_).receiver_noise_figure = SYS_config.UE_receiver_noise_figure;
    UEs(u_).BS_receiver_noise_figure = SYS_config.BS_receiver_noise_figure;
    
    % UL power control error,在这里产生是因为方便知道每个UE的power control error情况
    UEs(u_).PCe = normrnd(0,SYS_config.PCe_param);
    %     UEs(u_).PCe = -9+18*rand(1,1);
    
    % 接收功率热噪声（W）
    UEs(u_).thermal_noise_W_RB = 10^(0.1*SYS_config.UE.thermal_noise_density)/1000 * SYS_config.RB_bandwidth * 10^(UEs(u_).receiver_noise_figure/10);
    UEs(u_).BS_thermal_noise_W_RB = 10^(0.1*SYS_config.UE.thermal_noise_density)/1000 * SYS_config.RB_bandwidth * 10^(UEs(u_).BS_receiver_noise_figure/10);
    
    % 传入时钟模块
    UEs(u_).clock = networkClock;
end

%% 确保小区内有UE
if length(UEs)==1 && isempty(UEs(1).id)
    no_UEs = true;
else
    no_UEs = false;
end
BS_pos  = zeros(length(eNodeB_sites), 2);
for i = 1:length(eNodeB_sites)
    BS_pos(i, :) = eNodeB_sites(i).pos;
end

if ~no_UEs
   % 根据前面求得的基站服务的区域和UE的位置，给UE分配服务基站
    UE_positions_m = zeros(length(UEs),2);
    
    num_hetnet_sectors = networkPathlossMap.num_hetnet_sectors;
    if SYS_config.isDouble
        switch SYS_config.scene_type
            case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                num_first_UEs = (length(eNodeB_sectors)-num_hetnet_sectors)*SYS_config.UE_per_eNodeB;
            otherwise
                num_first_UEs = length(UEs)/2;
        end
    else
        num_first_UEs = 0;
    end
    networkPathlossMap.num_first_UEs = num_first_UEs;
    num_first_sectors = networkPathlossMap.num_first_sectors;
    for u_ = 1:length(UEs)     
        a = 1;
        % 根据基站服务的区域分配服务基站
        eNodeB_id = NR_calculate_attached_sector(SYS_config, BS_pos, UEs(u_).pos);
%         if theta(u_) <= 120
%            eNodeB_id = 31;
%         elseif theta(u_) <= 240
%            eNodeB_id = 32;
%         else
%            eNodeB_id = 33;
%         end
%         [ ~,~,eNodeB_id]    = networkPathlossMap.cell_assignment(SYS_config,u_,UEs(u_).pos,num_first_UEs);%将UE所在的扇区和eNB找出来
        
       % 将UE归附于相应的小区
        eNodeB_sectors(eNodeB_id).attachUser(UEs(u_));
        UE_positions_m(u_,:) = UEs(u_).pos;
        
        % 确定统计哪些UE
        if SYS_config.isWraparound
            compute_only_UEs_from_this_eNodeBs = 1:57;% 只统计被1到19号基站（1到57号小区）服务的UE，下面类似
        else
            if SYS_config.isDouble
                if SYS_config.isS2F % 异构中大场景到小场景与小场景到大场景
                    compute_only_UEs_from_this_eNodeBs = 1:num_first_sectors;% 求大场景对小场景的干扰
                else
                    compute_only_UEs_from_this_eNodeBs = num_first_sectors+1:length(eNodeB_sectors);% 求小场景对大场景的干扰
                end
            else
                compute_only_UEs_from_this_eNodeBs = 1:length(eNodeB_sectors);
            end
        end
        if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id == compute_only_UEs_from_this_eNodeBs,1))
            % 不统计该UE
            UEs(u_).deactivate_UE = true;%关闭UE统计
        else
            % 统计该UE
            UEs(u_).deactivate_UE = false;
        end
    end

end

function NR_init_RB_grid(SYS_config,eNodeB_sites)
% 初始化RB

% 为每个小区加入RB
for b_ = 1:length(eNodeB_sites)
    for s_=1:length(eNodeB_sites(b_).sectors)
        
        % 设置eNodeB是否为常开模式，无论是否有服务用户
        eNodeB_sites(b_).sectors(s_).always_on = SYS_config.always_on;
        
        % RB 分配表的创建与初始化
        eNodeB_sites(b_).sectors(s_).RB_grid = network_elements.resourceBlockGrid(SYS_config.N_RB);
        eNodeB_sites(b_).sectors(s_).RB_grid.set_homogeneous_power_allocation(SYS_config,eNodeB_sites(b_).sectors(s_).max_power,eNodeB_sites(b_).sectors(s_).signaling_power);
    end
end
