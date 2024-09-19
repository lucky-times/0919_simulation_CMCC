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

% %以中心点为基准，创建半径为SYS_config.UE_r，间隔为step度的UE飞行圆
% r = SYS_config.UE_r;
% step = 0.5;
% theta = 0:step:360-step;
% UE_positions = [r*cosd(theta); r*sind(theta)]';

str = '..\..\外场测试数据\7.4.2多层窄波束配置--蜂窝三扇区组网\7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度200米-加扰50%.csv';
data = readtable(str);
info = data(:, [3,4,5,8]);% 高度、经度、纬度、归属PCI
info = table2array(info);
temp = info((info(:, 1)>=190 & info(:, 1)<=215), :);% 高度、经度、纬度、归属PCI
r_100 = 92:148;
r_300 = 206:339;
r_500 = 377:604;
r_800 = 680:1045;
r_1200 = 1251:1804;
switch(SYS_config.UE_r)
    case 100
    r = r_100;
    case 300
    r = r_300;
    case 500
    r = r_500;
    case 800
    r = r_800;
    case 1200
    r = r_1200;
end




UE_height = temp(r, 1);
% 读取经度、纬度
longlat = temp(r, [2, 3]);
UE_positions = zeros(size(longlat));
figure;
hold on;
for i = 1:length(longlat)
    UE_positions(i, :) = latlon_to_xy(longlat(i, :));
    scatter(UE_positions(i, 1), UE_positions(i, 2));
end
axis equal;
title('UAV位置分布')

figure
plot(1:length(temp(r, 4)), temp(r, 4), LineWidth=2);
ylim([10, 310]);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(SYS_config.labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('PCI');
    str = sprintf('UAV Serving Cell (Field Data), h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
elseif(SYS_config.labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('PCI');
    str = sprintf('UAV服务小区（实测数据）, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
end
title(str)


%% 创建用户
UEs = network_elements.UE;
UE_pos_vector = zeros(size(UE_positions,1),2);
for u_ = 1:size(UE_positions,1)
    % 用户基本配置
    UEs(u_)     = network_elements.UE; %实体建立
    UEs(u_).id  = u_; %用户id设置
%     UEs(u_).height = SYS_config.UE_height;
    UEs(u_).height = UE_height(i);
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
    last_idx = [];
    for u_ = 1:length(UEs)     
        a = 1;
        % 根据基站服务的区域分配服务基站
        id_x = NR_calculate_attached_sector(SYS_config, eNodeB_sites,eNodeB_sectors, UEs(u_),last_idx);   
        last_idx = id_x;
       % 将UE归附于相应的小区
        eNodeB_sectors(id_x).attachUser(UEs(u_));
        UE_positions_m(u_,:) = UEs(u_).pos;
        
        % 确定统计哪些UE
%         if SYS_config.isWraparound
%             compute_only_UEs_from_this_eNodeBs = 1:57;% 只统计被1到19号基站（1到57号小区）服务的UE，下面类似
%         else
%             if SYS_config.isDouble
%                 if SYS_config.isS2F % 异构中大场景到小场景与小场景到大场景
%                     compute_only_UEs_from_this_eNodeBs = 1:num_first_sectors;% 求大场景对小场景的干扰
%                 else
%                     compute_only_UEs_from_this_eNodeBs = num_first_sectors+1:length(eNodeB_sectors);% 求小场景对大场景的干扰
%                 end
%             else
%                 compute_only_UEs_from_this_eNodeBs = 1:length(eNodeB_sectors);
%             end
%         end
%         if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id == compute_only_UEs_from_this_eNodeBs,1))
%             % 不统计该UE
%             UEs(u_).deactivate_UE = true;%关闭UE统计
%         else
%             % 统计该UE
%             UEs(u_).deactivate_UE = false;
%         end、
        UEs(u_).deactivate_UE = false;
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
