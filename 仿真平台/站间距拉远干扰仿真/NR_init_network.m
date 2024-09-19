function [eNodeB_sites,eNodeB_sectors,networkMacroscopicPathlossMap] = NR_init_network(SYS_config)
% 网络拓扑产生文件
% 输入参数：
% SYS_config：参数配置系统
% 输出参数：
% eNodeBs_sites：基站
% eNodeB_sectors：扇区
% networkMacroscopicPathlossMap：路损映射图
%

%% 部分数据赋值
data_res               = SYS_config.map_resolution;           % 地图分辨率
eNodeB_sector_tx_power = SYS_config.eNodeB_tx_power;          % 基站发射功率

%% 创建eNodeB实体
if SYS_config.debug_level>=1
    fprintf('Creating eNodeBs\n');
end
% 确定基站的位置
[eNodeB_sites,num_hetnet_sites] = NR_init_create_eNodeBs(SYS_config);
% 改变femto基站的位置，根据房间中心求房间区域
if SYS_config.isFemto
    for b_ = 1:length(eNodeB_sites)
        if strcmp(eNodeB_sites(b_).site_type,'indoor')
            room_centre = eNodeB_sites(b_).pos;% 原来的位置是房间的中心
            room_x = [room_centre(1)-10 room_centre(1)+10];
            room_y = [room_centre(2)-10 room_centre(2)+10];
            eNodeB_sites(b_).femto_room = [room_x;room_y];% 保存房间区域，画图的时候用
            new_pos = [random('unif',room_x(1),room_x(2)),random('unif',room_y(1),room_y(2))];% 给femto基站一个随机的位置
            eNodeB_sites(b_).pos = new_pos;
        end
    end
end

% 给基站其他参数赋值
s_idx   = 1;
eNodeB_sectors = network_elements.eNodeB_sector;  % 初始化,eNodeBs为eNodeB_sector的对象

for b_ = 1:length(eNodeB_sites)-num_hetnet_sites
   % 创建小区实体
    eNodeB_sites(b_).sectors    = network_elements.eNodeB_sector;
    for s_ = 1:length(SYS_config.sector_azimuths)
        eNodeB_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeB_sites(b_).sectors(s_).parent_eNodeB = eNodeB_sites(b_);
        eNodeB_sites(b_).sectors(s_).id            = s_; %配置小区id（相对于本eNodeB而言，取值为1,2,3）
        switch SYS_config.scene_type
            case {'UMA','RMa'}
                eNodeB_sites(b_).sectors(s_).azimuth       = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + SYS_config.sector_azimuths(s_));
            case {'UMI','UMa_to_UMi'}
                if ~SYS_config.isManhattan %Rand drop 的天线指向覆盖圈圆心
                    eNodeB_sites(b_).sectors(s_).azimuth = atan2(eNodeB_sites(b_).parent_centre_pos(2) -eNodeB_sites(b_).pos(2),eNodeB_sites(b_).parent_centre_pos(1)-eNodeB_sites(b_).pos(1))./pi*180;
                else
                    eNodeB_sites(b_).sectors(s_).azimuth = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + 330);
                end
            case {'InH','UMa_to_InH','UMi_to_InH','InH2','InH3'} %InH的天线偏角为0
                eNodeB_sites(b_).sectors(s_).azimuth = 0;
        end
        eNodeB_sites(b_).sectors(s_).max_power     = eNodeB_sector_tx_power; %配置基站的最大发射功率
        eNodeB_sites(b_).sectors(s_).antenna_type  = SYS_config.antenna.antenna_gain_pattern; % 天线类型
        
        eNodeB_sites(b_).sectors(s_).eNodeB_id     = s_idx; %设置小区的全局id号
        eNodeB_sectors(s_idx) = eNodeB_sites(b_).sectors(s_); %将上面是通过eNodeB_sites(b_).sectors(s_)进操作的，现在将其赋给eNodeB_sectors，以后可以通过eNodeB_sectors进行对扇区的操作
        
         % 为基站配置天线
        antennas.antenna.attach_antenna_to_eNodeB(eNodeB_sites(b_).sectors(s_),SYS_config,eNodeB_sites(b_).site_type);% 初始化天线
        
        % 创建宏路损模型
        if SYS_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);%b_为基站，s_为扇区
        end
        
        %双系统
        if SYS_config.isDouble 
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                    % 对异构来说，这步就是第一系统的
                    eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                        SYS_config,...
                        SYS_config.macroscopic_pathloss_model,...
                        SYS_config.frequency,...
                        SYS_config.macroscopic_pathloss_model_settings);
                otherwise
                    if b_<=length(eNodeB_sites)/2 %第一个系统路损模型
                        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                            SYS_config,...
                            SYS_config.macroscopic_pathloss_model,...
                            SYS_config.frequency,...
                            SYS_config.macroscopic_pathloss_model_settings);
                    else %第二个系统路损模型
                        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                            SYS_config,...
                            SYS_config.macroscopic_pathloss_model2,...
                            SYS_config.frequency2,...
                            SYS_config.macroscopic_pathloss_model_settings2);
                    end
            end
        else %单系统的路损模型
            eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                SYS_config,...
                SYS_config.macroscopic_pathloss_model,...
                SYS_config.frequency,...
                SYS_config.macroscopic_pathloss_model_settings);
        end
        s_idx = s_idx + 1;
    end
end
num_hetnet_sectors = 0;% 异构场景中，第二个系统的扇区数
for b_ = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites)% 若不是异构网络这段不执行,num_hetnet_sites~=0的话表明是异构网络
    eNodeB_sites(b_).sectors = network_elements.eNodeB_sector;
    for s_ = 1:length(SYS_config.sector_azimuths2)
        eNodeB_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeB_sites(b_).sectors(s_).parent_eNodeB = eNodeB_sites(b_);
        eNodeB_sites(b_).sectors(s_).id            = s_;      
        switch eNodeB_sites(b_).site_type
            case 'macro'
                eNodeB_sites(b_).sectors(s_).azimuth = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + SYS_config.sector_azimuths2(s_));
            case 'micro'
                eNodeB_sites(b_).sectors(s_).azimuth = atan2(eNodeB_sites(b_).parent_centre_pos(2) -eNodeB_sites(b_).pos(2),eNodeB_sites(b_).parent_centre_pos(1)-eNodeB_sites(b_).pos(1))./pi*180;
        end
        eNodeB_sites(b_).sectors(s_).max_power     = SYS_config.eNodeB_tx_power2;
        eNodeB_sites(b_).sectors(s_).antenna_type  = SYS_config.antenna.antenna_gain_pattern;
        eNodeB_sites(b_).sectors(s_).eNodeB_id     = s_idx;
        eNodeB_sectors(s_idx)                              = eNodeB_sites(b_).sectors(s_);
        antennas.antenna.attach_antenna_to_eNodeB(eNodeB_sites(b_).sectors(s_),SYS_config,eNodeB_sites(b_).site_type);
        
        if SYS_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);%b_为基站，s_为扇区
        end
        
        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            SYS_config,...
            SYS_config.macroscopic_pathloss_model2,...
            SYS_config.frequency2,...
            SYS_config.macroscopic_pathloss_model_settings2);
        
        s_idx = s_idx + 1;
        num_hetnet_sectors = num_hetnet_sectors + 1;
    end
end

%% 确定地图大小
switch SYS_config.scene_type
    case {'UMA','RMa','UMa_to_InH'}
        % 得到基站的位置
        tx_pos = zeros(length(eNodeB_sites),2);%2列，第一列为x，第二列为y
        for b_ = 1:length(eNodeB_sites)
            tx_pos(b_,:) = eNodeB_sites(b_).pos;
        end
        
        % 根据基站的位置计算地图边界
        if  SYS_config.nr_eNodeB_rings ~= 0
            roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
            roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
        else % 如果基站只有一个
            roi_x = [-SYS_config.ISD,SYS_config.ISD];
            roi_y = [-SYS_config.ISD,SYS_config.ISD];
        end
        
        ISD = SYS_config.ISD;
        roi_x =[roi_x(1)-2*ISD/3 roi_x(2)+2*ISD/3];
        roi_y =[roi_y(1)-2*ISD/3 roi_y(2)+2*ISD/3];
    case 'UMI'
        if ~SYS_config.isManhattan
            % 根据ISD算出来的，能刚好包括全部基站及基站的覆盖范围
            ISD = SYS_config.ISD;
            roi_x =[-(ISD*2)-(ISD*2)/3 (ISD*2)+(ISD*2)/3];
            roi_y =[-(ISD*(3^0.5))-(ISD*2)/3 (ISD*(3^0.5))+(ISD*2)/3];
        else
            % 固定的大小
            roi_x =[0 1080];
            roi_y =[0 1080];
        end
    case 'UMa_to_UMi'
        if ~SYS_config.isManhattan
            ISD = SYS_config.ISD;
            roi_x =[-(ISD*2)-(ISD*2)/3 (ISD*2)+(ISD*2)/3];
            roi_y =[-(ISD*(3^0.5))-(ISD*2)/3 (ISD*(3^0.5))+(ISD*2)/3];
        else
            tx_pos = zeros(length(eNodeB_sites),2);
            for b_ = 1:length(eNodeB_sites)
                tx_pos(b_,:) = eNodeB_sites(b_).pos;
            end
            roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
            roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
            ISD = SYS_config.ISD;
            roi_x =[roi_x(1)-2*ISD/3 roi_x(2)+2*ISD/3];
            roi_y =[roi_y(1)-2*ISD/3 roi_y(2)+2*ISD/3];
        end
    case 'UMi_to_InH'
        UMi_r = SYS_config.UMi_r;
        tx_pos = zeros(length(eNodeB_sites),2);
        for b_ = 1:length(eNodeB_sites)
            tx_pos(b_,:) = eNodeB_sites(b_).pos;
        end
        roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
        roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
        roi_x =[roi_x(1)-UMi_r roi_x(2)+UMi_r];
        roi_y =[roi_y(1)-UMi_r roi_y(2)+UMi_r];
    case {'InH','InH2'}
        % 固定的大小
        roi_x =[0 120];
        if ~SYS_config.isFemto
            roi_y =[0 50];
        else
            roi_y =[0 40];
        end
    case {'InH3'}
        % 固定的大小
        roi_x =[0 120];
        if ~SYS_config.isFemto
            roi_y =[0 100+SYS_config.InH3_d];
        else
            roi_y =[0 40];
        end
end

%% 初始化路损图谱
networkMacroscopicPathlossMap                        = channel_gain_wrappers.macroscopicPathlossMap;%无参构造函数
networkMacroscopicPathlossMap.data_res               = data_res;
networkMacroscopicPathlossMap.roi_x                  = roi_x;
networkMacroscopicPathlossMap.roi_y                  = roi_y;
networkMacroscopicPathlossMap.num_hetnet_sites       = num_hetnet_sites;
networkMacroscopicPathlossMap.num_hetnet_sectors     = num_hetnet_sectors;

if SYS_config.isDouble
    switch SYS_config.scene_type
        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
            num_first_sites = length(eNodeB_sites) - num_hetnet_sites;
            num_first_sectors = length(eNodeB_sectors)-num_hetnet_sectors;
        otherwise
            num_first_sites = length(eNodeB_sites)/2;
            num_first_sectors = length(eNodeB_sectors)/2;
    end
else
    % 单系统通过诸如length（eNodeB_sites）调用较快
    num_first_sites = 0;
    num_first_sectors = 0;
end
networkMacroscopicPathlossMap.num_first_sectors = num_first_sectors;% 双系统中第一个系统的扇区数
networkMacroscopicPathlossMap.num_first_sites = num_first_sites;% 双系统中第一个系统的基站数


%% 计算路径损耗
if SYS_config.debug_level>=1
    fprintf('Creating cell pathloss map\n');
end
% macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(SYS_config,eNodeB_sites,networkMacroscopicPathlossMap);

%% 从以前的平台继承来的，没什么用
networkMacroscopicPathlossMap.name = 'NR or LTE pathloss model';
