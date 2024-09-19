classdef generalPathlossModel < handle
% 产生路损模型与天线单元增益
    properties
        name
    end
    methods (Abstract)
        pathloss_in_db = pathloss(distance)
    end
    methods (Static)
        function calculate_pathloss_maps(SYS_config,eNodeB_sites,networkMacroscopicPathlossMap)
            % 给一些变量初始化
            total_sectors         = length([eNodeB_sites.sectors]);
            data_res              = networkMacroscopicPathlossMap.data_res;
            roi_x                 = networkMacroscopicPathlossMap.roi_x;
            roi_y                 = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels    = NR_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels     = roi_maximum_pixels(2);
            roi_width_pixels      = roi_maximum_pixels(1);
            distance_matrix       = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
            cell_pathloss_data    = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            sector_antenna_gain   = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            % 与三维地图有关变量
            UE_antenna_gain_1     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            UE_antenna_gain_2     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            UE_antenna_gain     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            penetration_loss = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            d_2d_in = zeros(roi_height_pixels,roi_width_pixels);
            
            site_positions     = reshape([eNodeB_sites.pos],2,[])';
            site_positions_pix = NR_common_pos_to_pixel( site_positions, [roi_x(1) roi_y(1)], data_res);
            %% 创建3维地图
            switch SYS_config.scene_type
                case {'UMA' , 'UMI','UMa_to_UMi'}
                    [row,col] = size(d_2d_in);
                    if ~SYS_config.isManhattan %UMA及UMi的Rand drop模型
                        % 'UMI','UMa_to_UMi'中微小区不是曼哈顿执行这一段
%                         if ~SYS_config.isNewUMA %先前的UMA场景（500ISD）的室内外比例为20%、20%
%                             index = NR_randperm(row*col,round(row*col*0.2));% 20%在室内
%                         else
%                             index = NR_randperm(row*col,round(row*col*0));% 0%在室内
%                         end
                        switch SYS_config.scene_type
                            case 'UMA'
                                index = NR_randperm(row*col,round(row*col*0));
                            case 'UMI'
                                index = NR_randperm(row*col,round(row*col*0.2));
                            case 'UMa_to_UMi'
                                index = NR_randperm(row*col,round(row*col*0.2));
                        end
                        for i = 1:length(index)
                            if index == 0
                                d_2d_in = zeros(size(d_2d_in));
                            else
                                d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1)); %根据d_2d_in（距离室外的距离）区分室内外用户，并进行路损计算
                            end
                        end
                    else % 曼哈顿
                        % 'UMI','UMa_to_UMi'中微小区是曼哈顿执行这一段
                        switch SYS_config.scene_type
                            case 'UMI'
                                % 'UMI'中微小区是曼哈顿
                                % 若地图分辨率为5，则一个像素点占的距离为5m
                                % 因为楼宽和街道宽的特殊性，SYS_config.map_resolution应该是一个能整除75和15的数
                                % 若不能整除，对于曼哈顿系统自己的共存没有影响，但是对宏蜂窝与曼哈顿的异构会报错，因为矩阵规模不一致
                                n_ = fix(75/SYS_config.map_resolution);% 楼宽
                                m_ = fix(15/SYS_config.map_resolution);% 街道宽
                                % 从左下角开始往上添加建筑物
                                for row_=m_+1:m_+n_ % 第一行建筑物
                                    d_2d_in(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                    for i=1:11
                                        d_2d_in(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));%产生两个随机值并选择较小的
                                    end
                                end
                                for j=1:11% 其余11行建筑物
                                    for row_=m_+j*n_+j*m_+1:m_+(j+1)*n_+j*m_
                                        d_2d_in(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                        for i=1:11
                                            d_2d_in(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));
                                        end
                                    end
                                end
                            case 'UMa_to_UMi'
                                % 'UMa_to_UMi'中微小区是曼哈顿
                                % UMa场景下的d_2d_in
                                if ~SYS_config.isNewUMA %先前的UMA场景（500ISD）的室内外比例为80%、20%
                                    index = NR_randperm(row*col,round(row*col*0.8));% 80%在室内
                                else
                                    index = NR_randperm(row*col,round(row*col*0.2));% 20%在室内
                                end
                                for i = 1:length(index)
                                    d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1)); %根据d_2d_in（距离室外的距离）区分室内外用户，并进行路损计算
                                end
                                % 曼哈顿场景下的d_2d_in
                                O_M = [-540,-540];
                                K_M = [540,540];
                                O_M_pix = NR_common_pos_to_pixel(O_M, [roi_x(1) roi_y(1)], data_res);
                                K_M_pix = NR_common_pos_to_pixel(K_M, [roi_x(1) roi_y(1)], data_res);
                                M_row = O_M_pix(2):K_M_pix(2);
                                M_col = O_M_pix(1):K_M_pix(1);
                                d_2d_in_M = zeros(size(d_2d_in(M_row,M_col)));
                                
                                n_ = round(75/SYS_config.map_resolution);
                                m_ = round(15/SYS_config.map_resolution);
                                % 根据曼哈顿拓扑，如下8个*号组成一个建筑物
                                % 添加建筑物的方法就是一行一行像素点地添加
                                % 一行建筑物的意思：如下图，就是由两行像素点（*）构成一行建筑物
                                % ****  ****
                                % ****  ****
                                % 
                                % ****  ****
                                % ****  ****
                                % 从左下角开始往上添加建筑物
                                for row_=m_+1:m_+n_ % 第一行建筑物
                                    d_2d_in_M(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                    for i=1:11
                                        d_2d_in_M(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));%产生两个随机值并选择较小的
                                    end
                                end
                                for j=1:11% 其余11行建筑物
                                    for row_=m_+j*n_+j*m_+1:m_+(j+1)*n_+j*m_
                                        d_2d_in_M(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                        for i=1:11
                                            d_2d_in_M(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));
                                        end
                                    end
                                end
                                % 先求出UMa场景下的d_2d_in，然后用曼哈顿下的d_2d_in_M替换d_2d_in对应的部分
                                d_2d_in(M_row,M_col) = d_2d_in_M;
                        end
                    end
                    % ue高度
                    h_ut = ones(size(d_2d_in)) * SYS_config.UE_height;
                    
                    % 穿损——一般低穿损，一半低穿损

                    % 低穿损 （根据38.900里的穿损公式计算）
                    L_glass = 2+0.2*SYS_config.frequency/1e9;
                    L_concrete = 5+4*SYS_config.frequency/1e9;
                    N_low = 0+4.4*randn(row,col);
                    low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                    index = NR_randperm(row*col,round(row*col*0.5));%低路损、高路损各取一半
                    low(index) = 0;
                    % 高穿损
                    L_IRRglass = 23+0.3*SYS_config.frequency/1e9;
                    N_high = 0+6.5*randn(row,col);
                    high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*d_2d_in+N_high;
                    [m,n]=find(low(:,:)~=0);
                    for i = 1:length(m)
                        high(m(i),n(i)) = 0;
                    end
                    penetration = low+high;%总的传损
                    [m,n]=find(d_2d_in(:,:)==0);
                    for i = 1:length(m)
                        penetration(m(i),n(i)) = 0; %室外用户穿损为0
                    end
                    if SYS_config.isDouble% 频率不一样，穿损不一样
                        % 低穿损
                        L_glass = 2+0.2*SYS_config.frequency2/1e9;
                        L_concrete = 5+4*SYS_config.frequency2/1e9;
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                        low(index) = 0;
                        %高穿损
                        L_IRRglass = 23+0.3*SYS_config.frequency2/1e9;
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*d_2d_in+N_high;
                        [m,n]=find(low(:,:)~=0);
                        for i = 1:length(m)
                            high(m(i),n(i)) = 0;
                        end
                        penetration2 = low+high;
                        [m,n]=find(d_2d_in(:,:)==0);
                        for i = 1:length(m)
                            penetration2(m(i),n(i)) = 0;
                        end
                    end
                case 'RMa'
                    % 50%室内，50%车上，室内低穿损，车上没有穿损只有一个N(μ, σP）
                    [row,col] = size(d_2d_in);
                    index = NR_randperm(row*col,round(row*col*0.5));% 50%在室内
                    for i = 1:length(index)
                        d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1));
                    end
                    % ue高度
                    h_ut = ones(size(d_2d_in)) * 1.5;
                    
                    % 低穿损                    
                    L_glass = 2+0.2*SYS_config.frequency/1e9;
                    L_concrete = 5+4*SYS_config.frequency/1e9;
                    N_low = 0+4.4*randn(row,col);
                    low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                    penetration = low;
                    [m,n]=find(d_2d_in(:,:)==0);
                    for i = 1:length(m)
                        penetration(m(i),n(i)) = 9+5*randn;% 车上没有穿损只有一个N(μ, σP）
                    end
                    if SYS_config.isDouble
                        L_glass = 2+0.2*SYS_config.frequency2/1e9;
                        L_concrete = 5+4*SYS_config.frequency2/1e9;
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                        penetration2 = low;
                        [m,n]=find(d_2d_in(:,:)==0);
                        for i = 1:length(m)
                            penetration2(m(i),n(i)) = 9+5*randn;
                        end
                    end
            end
            
            % 初始化距离和角度矩阵
            position_grid_pixels      = zeros(roi_height_pixels*roi_width_pixels,2);
            position_grid_pixels(:,1) = reshape(repmat(1:roi_width_pixels,roi_height_pixels,1),1,roi_width_pixels*roi_height_pixels);
            position_grid_pixels(:,2) = repmat(1:roi_height_pixels,1,roi_width_pixels)';
            position_grid_meters      = NR_common_pixel_to_pos(position_grid_pixels,networkMacroscopicPathlossMap.coordinate_origin,networkMacroscopicPathlossMap.data_res);
            
            % 求femto的qcol和qrow，qcol是经过垂直的墙的个数，qrow是经过水平的墙的个数
            if SYS_config.isFemto
                if ~isFemto_het
                    q_f_matrix = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
                    % q_col
                    % 每个房间的大小是20*20，于是可以通过这种算法判断该点在哪个房间
                    q_col_mod = mod(position_grid_meters(:,1),20);
                    index = find(q_col_mod == 0);
                    q_col = fix(position_grid_meters(:,1)/20)+1;
                    q_col(index) = q_col(index) - 1;
                    % q_row
                    q_row_mod = mod(position_grid_meters(:,2),20);
                    index = find(q_row_mod == 0);
                    q_row = fix(position_grid_meters(:,2)/20)+1;
                    q_row(index) = q_row(index) - 1;
                else
                    q_f_matrix = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
                    % q_col
                    position_col_matrix = reshape(position_grid_meters(:,1),roi_height_pixels,roi_width_pixels);
                    insite_col = position_col_matrix(InH_row,InH_col) + 60;% 加60是为了调整坐标，使其和非异构的坐标一致，方便下面计算
                    q_col_mod = mod(insite_col,20);
                    index = find(q_col_mod == 0);
                    q_col = fix(insite_col/20)+1;
                    q_col(index) = q_col(index) - 1;
                    % q_row
                    position_row_matrix = reshape(position_grid_meters(:,2),roi_height_pixels,roi_width_pixels);
                    insite_row = position_row_matrix(InH_row,InH_col) + 20;
                    q_row_mod = mod(insite_row,20);
                    index = find(q_row_mod == 0);
                    q_row = fix(insite_row/20)+1;
                    q_row(index) = q_row(index) - 1;
                end
            end
            
            %% 计算小区的路损
            s_idx = 1;
            all_sectors   = [eNodeB_sites.sectors];
            eNodeB_id_set = [all_sectors.eNodeB_id];
            % 2panel,+-180,UE随机方向角
            rnd_phi_panel1=180-360*rand(size(h_ut)); %加上-180~180的随机旋转相位
            rnd_phi_panel2=rnd_phi_panel1-180;%两个panel差180度
            for b_ = 1:length(eNodeB_sites)
                distances               = sqrt((position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)).^2 + (position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)).^2);%每个像素点到enb的距离
                % 按地图的size调整距离矩阵
                distance_matrix(:,:,b_) = reshape(distances,roi_height_pixels,roi_width_pixels);
                
                % 求像素点距离基站隔了多少面墙
                if SYS_config.isFemto
                    % 先求基站在哪个房间
                    if ~isFemto_het
                        site_pos = eNodeB_sites(b_).pos;
                    else
                        site_pos = [eNodeB_sites(b_).pos(1)+60,eNodeB_sites(b_).pos(2)+20];% 异构场景下需要调整坐标使其和非异构的坐标一致
                    end
                    if mod(site_pos(1),20) == 0
                        site_col =  site_pos(1)/20;
                    else
                        site_col = fix(site_pos(1)/20)+1;
                    end
                    if mod(site_pos(2),20) == 0
                        site_row =  site_pos(2)/20;
                    else
                        site_row = fix(site_pos(2)/20)+1;
                    end
                    col_f = abs(site_col - q_col);
                    row_f = abs(site_row - q_row);
                    q_f = col_f + row_f;% 垂直墙和水平墙加起来为总的墙
                    if ~isFemto_het
                        q_f_matrix(:,:,b_) = reshape(q_f,roi_height_pixels,roi_width_pixels);
                    else
                        if strcmp(eNodeB_sites(b_).site_type, 'indoor')
                            % 只有indoor基站才有这个值，其他类型的基站为0
                            q_f_matrix(InH_row,InH_col,b_) = q_f;
                        end
                    end
                end
                
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    
                    % 避免像素点到基站的距离小于室内距离的情况,这时室内距离等于像素点到基站距离的一半
                    d_2d_tmp = distance_matrix(:,:,b_);
                    d_2d_in(find(d_2d_tmp(:,:)<d_2d_in(:,:))) = d_2d_tmp(find(d_2d_tmp(:,:)<d_2d_in(:,:)))/2;
                    
                    %调用pathloss计算模块，得到路径损耗与LOS分布
                    [cell_pathloss_data(:,:,s_idx),LOS_compare(:,:,s_idx)] = eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model.pathloss(distance_matrix(:,:,b_),h_ut,d_2d_in);
                    
                    % 将穿透损耗表示为三维矩阵的形式方便先pathloss一样调用
                    if ~SYS_config.isDouble
                        penetration_loss(:,:,s_idx) = penetration;
                    else
                        num_first_sectors = networkMacroscopicPathlossMap.num_first_sectors;
                        if s_idx <= num_first_sectors
                            penetration_loss(:,:,s_idx) = penetration;
                        else
                            penetration_loss(:,:,s_idx) = penetration2;
                        end
                    end
                    if SYS_config.isFemto
                        % femto中每隔一面墙穿透损耗5dB
                        % 对非异构来说，右边的penetration_loss为0
                        % 对异构来说，右边的penetration_loss(:,:,第二系统基站)表示室外基站到室内的穿墙损耗，q_f_matrix(:,:,第二系统基站)=0
                        penetration_loss(:,:,s_idx) = penetration_loss(:,:,s_idx) + q_f_matrix(:,:,b_) * 5;
                    end
                    
                    %传入基站高度参数和最大天线增益
                    num_first_sites = networkMacroscopicPathlossMap.num_first_sites;
                    switch SYS_config.scene_type
                        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                            if b_ <= num_first_sites
                                eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height;
                                eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain;
                            else
                                eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height2;
                                eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain2;
                            end
                        otherwise
                            eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height;
                            eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain;
                    end
                    % 计算水平偏角
                    switch eNodeB_sites(b_).site_type
                        case 'indoor'
                            angle_grid = (180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),eNodeB_sites(b_).sectors(s_).tx_height-SYS_config.UE_height)); %等效坐标系旋转90度，左右方向（x）表示theta，前后方向（y）表示phi。这里phi以垂直向下方向为参考系
                            angle_grid_ue=180-(180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)));
                        otherwise
                            angle_grid = (180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),(position_grid_meters(:,1)-eNodeB_sites(b_).pos(1))))-eNodeB_sites(b_).sectors(s_).azimuth;
                            angle_grid_ue=180-angle_grid;%UE和BS的水平角互补
                    end
                    horizontal_angle_grid_ue   = reshape(angle_grid_ue,roi_height_pixels,roi_width_pixels);
                    horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(horizontal_angle_grid_ue);
                    %转换为(-180,180]
                    horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
                    if strcmp(SYS_config.antenna.antenna_gain_pattern,'NRAntennaBeamforming')
                        
                        % 水平偏转角等效到-180~180
                        horizontal_angle_grid   = reshape(angle_grid,roi_height_pixels,roi_width_pixels);
                        horizontal_angle_grid_s = utils.miscUtils.wrapTo359(horizontal_angle_grid);
                        horizontal_angle_grid_s = utils.miscUtils.wrapToAll180(horizontal_angle_grid_s);
                        
                        % 垂直偏角角度计算
                        switch eNodeB_sites(b_).site_type
                            case 'indoor'
                                vertical_angle_grid_el_ue=(180/pi)*(atan2(distance_matrix(:,:,b_),eNodeB_sites(b_).sectors(s_).tx_height - SYS_config.UE_height));
                                vertical_angle_grid_el = (180/pi)*(atan2( sqrt((SYS_config.site_height - SYS_config.UE_height).^2+(position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)).^2),position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)));%theta以水平朝右为参考系
                                vertical_angle_grid_el   = reshape(vertical_angle_grid_el,roi_height_pixels,roi_width_pixels);
                            otherwise
                                vertical_angle_grid_el = (180/pi)*(atan2( distance_matrix(:,:,b_),h_ut(:,:)-eNodeB_sites(b_).sectors(s_).tx_height));
                                vertical_angle_grid_el_ue=(180/pi)*(atan2( distance_matrix(:,:,b_),eNodeB_sites(b_).sectors(s_).tx_height-h_ut(:,:)));
                        end
                        
                        phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
                        phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
                        %加上两侧panel随机旋转角后仍要保证总体方位角在-180到180间
                        phi_1 = utils.miscUtils.wrapToAll180(phi_1);
                        phi_2 = utils.miscUtils.wrapToAll180(phi_2);
                        UE_test = network_elements.UE;%用于调用element gain函数
                        antennas.antenna.attach_antenna_to_UE(UE_test,SYS_config);
                        
                        switch eNodeB_sites(b_).site_type
                            case 'macro'
                                %计算两个panel对应的水平旋转角，并选择天线增益最大的UE偏转角

                                sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                UE_antenna_gain(:,:,s_idx)  = max(UE_antenna_gain_1(:,:,s_idx),UE_antenna_gain_2(:,:,s_idx));
                                
                                switch SYS_config.attatch_mode%设置接入模式
                                    
                                    case 0 % 只通过路损接入
                                        sector_antenna_gain(:,:,s_idx)=0;
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 1 %通过路损与BS端的天线单元增益接入
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 2 %通过路损与UE端的天线单元增益接入
                                        sector_antenna_gain(:,:,s_idx)=0;
                                    case 3 %通过路损与BS与UE的天线单元增益接入
                                end
                            case 'micro'
                                if ~SYS_config.isManhattan
                                    sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s);
                                else
                                    sector_antenna_gain(:,:,s_idx) = 6;
                                end
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                UE_antenna_gain  = max(UE_antenna_gain_1,UE_antenna_gain_2);
                                
                                switch SYS_config.attatch_mode
                                    case 0
                                        sector_antenna_gain(:,:,s_idx)=0;
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 1
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 2
                                        sector_antenna_gain(:,:,s_idx)=0;
                                    case 3
                                end
                                
                            case 'indoor'
                                sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s);
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                
                                UE_antenna_gain  = max(UE_antenna_gain_1,UE_antenna_gain_2);
                                switch SYS_config.attatch_mode
                                    case 0
                                        sector_antenna_gain = zeros(size(sector_antenna_gain));
                                        UE_antenna_gain=zeros(size(UE_antenna_gain));
                                    case 1
                                        
                                        UE_antenna_gain=zeros(size(UE_antenna_gain));
                                    case 2
                                        sector_antenna_gain=zeros(size(sector_antenna_gain));
                                    case 3
                                end
                        end
                    else %如果使用其他天线模型
                        error('"%s" antenna is not supported',SYS_config.antenna.antenna_gain_pattern);
                    end
                    
                    networkMacroscopicPathlossMap.sector_idx_mapping(eNodeB_sites(b_).sectors(s_).eNodeB_id,:) = [b_ s_];
                    networkMacroscopicPathlossMap.site_sector_mapping(b_,s_)  = eNodeB_sites(b_).sectors(s_).eNodeB_id;
                    
                    s_idx = s_idx + 1;
                end
            end
            
            %% 微小区中心到各像素点距离
            micro_centre_distance_matrix = zeros(size(distance_matrix));
            for b_ = 1:length(eNodeB_sites)
                if strcmp(eNodeB_sites(b_).site_type,'micro')
                    if ~SYS_config.isManhattan
                        micro_centre__distances = sqrt((position_grid_meters(:,1)-eNodeB_sites(b_).parent_centre_pos(1)).^2 + (position_grid_meters(:,2)-eNodeB_sites(b_).parent_centre_pos(2)).^2);
                        micro_centre_distance_matrix(:,:,b_) = reshape(micro_centre__distances,roi_height_pixels,roi_width_pixels);% 微小区中心到地图上各个像素点的距离
                    end
                end
            end

            %% 保存数据
            cell_pathloss_data(isnan(cell_pathloss_data) | (cell_pathloss_data<0)) = 0;
            
            networkMacroscopicPathlossMap.micro_centre_distance_matrix = micro_centre_distance_matrix;% 微小区中心到地图上各个像素点的距离
            networkMacroscopicPathlossMap.distance_matrix = distance_matrix;% 基站到各像素点距离
            networkMacroscopicPathlossMap.H_height = h_ut;% 像素点高度
            networkMacroscopicPathlossMap.H_orientation = rnd_phi_panel1;% 旋转角
            networkMacroscopicPathlossMap.penetration_loss = penetration_loss;% 穿损，和路损是一样规模的矩阵
            networkMacroscopicPathlossMap. pathloss_900 = cell_pathloss_data;% 不包括天线增益的路损
            % 计算出的路损
            networkMacroscopicPathlossMap.pathloss(:,:,eNodeB_id_set) = cell_pathloss_data - sector_antenna_gain-UE_antenna_gain+SYS_config.cable_loss;% 包括天线增益的路损
            networkMacroscopicPathlossMap.LOS_compare(:,:,eNodeB_id_set) = LOS_compare;% LOS的概率矩阵
            networkMacroscopicPathlossMap.BS_antenna_gain = sector_antenna_gain;% 基站的天线增益
            networkMacroscopicPathlossMap.UE_antenna_gain = UE_antenna_gain;
        end
        
        function macroscopic_pathloss_model = generateMacroscopicPathlossModel(SYS_config,macroscopic_pathloss_model_name,frequency,macroscopic_pathloss_model_settings)
            % 通过macroscopic_pathloss_model_name参数采用合适的路损模型
            print_output = true;
            switch macroscopic_pathloss_model_name
                case 'TS36942'% LTE路损模型
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS36942PathlossModel(frequency,macroscopic_pathloss_model_settings.environment);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 36.942-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end
                case 'TS38900'% 5G路损模型
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,macroscopic_pathloss_model_settings.environment,SYS_config);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 38.900-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end
                case 'TS38901'% 5G路损模型
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,macroscopic_pathloss_model_settings.environment,SYS_config);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 38.901-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end 
                otherwise
                    error('"%s" macroscopic pathloss model not supported',macroscopic_pathloss_model_name);
            end
        end
    end
end
