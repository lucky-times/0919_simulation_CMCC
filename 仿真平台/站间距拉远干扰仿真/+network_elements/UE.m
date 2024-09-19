classdef UE < handle
    % UE实体文件
    % 具体功能：用于计算SINR

    properties

        id                    % UE的id
        pos                   % UE的位置(x,y)
        attached_site         % UE归属的eNode站址
        attached_sector_idx   % UE归属的eNode站址·
        attached_eNodeB       % UE归属的eNodeB
        height                % UE的高度
        orientation           % UE的旋转角
        antenna               % UE归属的eNodeB
        walking_model         % UE的运动模型
        downlink_channel      % UE的下行信道
        RB_grid               % UE的RB分配表
        ul_ad_TX_power   %UE上行功控后的发射功率

        uplink_channel        % UE的上行信道
        PCe                   % 上行功控误差
        receiver_noise_figure % UE的接收噪声系数
        BS_receiver_noise_figure % 上行时BS接收端的噪声
        thermal_noise_W_RB    % 热噪声（单位W）
        BS_thermal_noise_W_RB
        UE_tx_power           % UE上行时的发射功率
        penetration_loss      % 该UE到各个基站的穿透损耗
        rx_antenna_gain       % UE接收天线增益
        beam_antenna_gain     % UE波束赋型增益

        clock                 % 网络时钟

        wideband_SINR        % 保存下行计算的SINR
        wideband_SINR2       % 保存采用论文所提方法计算出SINR
        ul_wideband_SINR     % 保存上行计算的SINR
        ul_dl_wideband_SINR  % 保存上行干扰下行时的SINR
        dl_ul_wideband_SINR  % 保存下行干扰上行时的SINR
        bs_bs_pathloss       % 基站对基站的耦合损耗
        ue_ue_pathloss       % ue对ue的耦合损耗

        dl_ul_UE_gain        %下行干扰上行的被干扰链路的UE发射功率
        dl_ul_BS_gain        %下行干扰上行的被干扰链路的BS发射功率

        ul_dl_UE_gain        %上行干扰下行的被干扰链路的UE发射功率
        ul_dl_BS_gain        %下行干扰上行的被干扰链路的BS发射功率


        SINR_TTI             % 保存每个TTI中的下行SINR
        ul_SINR_TTI          % 保存每个TTI中的上行SINR
        ul_dl_SINR_TTI       % 保存每个TTI中的上行干扰下行SINR
        dl_ul_SINR_TTI       % 保存每个TTI中的下行干扰上行SINR
        attached_eNodeB_id_TTI % 保存每个TTI中UE连接的基站
        deactivate_UE        % 统计该UE与否
        cell_change          % 结构体与小区切换有关

        interfering_UE
        interfering_BS

    end

    methods
        % 初始化构造函数
        function obj = UE
            obj.deactivate_UE             = false;
            obj.cell_change.requested     = false;
            obj.cell_change.target_eNodeB = [];
        end

        function print(obj)
            if isempty(obj.attached_site)
                fprintf('User %d, (%d,%d), not attached to an eNodeB\n',obj.id,obj.pos(1),obj.pos(2));
            else
                fprintf('User %d, (%d,%d), Site %d, sector %d (eNodeB %d)\n',obj.id,obj.pos(1),obj.pos(2),obj.attached_site.id,obj.attached_sector_idx,obj.attached_eNodeB);
            end
            obj.walking_model.print;
        end

        %% 清空变量
        function clear(obj)
            obj.attached_site             = [];
            obj.attached_eNodeB           = [];
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.clock                     = [];
        end

        %% 根据移动模型进行移动
        function move(obj)
            new_pos = obj.walking_model.move(obj.pos);
            obj.pos = new_pos;
        end

        function move_back(obj)
            old_pos = obj.walking_model.move_back(obj.pos);
            obj.pos = old_pos;
        end
        %% 判断UE是否在ROI内
        function UE_in_roi = is_in_roi(obj,config,roi_x_range,roi_y_range)

            roi_x = obj.downlink_channel.macroscopic_pathloss_model.roi_x;
            roi_y = obj.downlink_channel.macroscopic_pathloss_model.roi_y;
            data_res = obj.downlink_channel.macroscopic_pathloss_model.data_res;
            sector_assignment = obj.downlink_channel.macroscopic_pathloss_model.sector_assignment;
            sector_assignment_double = obj.downlink_channel.macroscopic_pathloss_model.sector_assignment_double;
            num_first_UEs = obj.downlink_channel.macroscopic_pathloss_model.num_first_UEs;

            UE_pos_pix = NR_common_pos_to_pixel( obj.pos, [roi_x(1) roi_y(1)], data_res);

            try % 利用数组越界的错误判断UE代表的像素点是否在ROI中
                if config.isDouble
                    if UE_id<=num_first_UEs
                        s_idx = sector_assignment(UE_pos_pix(2),UE_pos_pix(1));
                    else
                        s_idx = sector_assignment_double(UE_pos_pix(2),UE_pos_pix(1));
                    end
                else
                    s_idx = sector_assignment(UE_pos_pix(2),UE_pos_pix(1));
                end
            catch
                UE_in_roi = false;
                return;
            end

            % 通过像素点不能完全判断UE的新位置是否在ROI内,如inh中(120.4750,27.1949)就判断在ROI内
            UE_pos_x = obj.pos(1);
            UE_pos_y = obj.pos(2);
            if UE_pos_x<roi_x_range(1) || UE_pos_x>roi_x_range(2)
                UE_in_roi = false;
                return;
            end
            if UE_pos_y<roi_y_range(1) || UE_pos_y>roi_y_range(2)
                UE_in_roi = false;
                return;
            end

            if s_idx == -1 % UE在ROI上但不在基站的服务范围内，在展示用的GUI上表现为白色的部分
                UE_in_roi = false;
                return;
            else
                UE_in_roi = true;
            end

            if UE_in_roi
                num_first_UEs = obj.downlink_channel.macroscopic_pathloss_model.num_first_UEs;
                [~,~,new_eNodeB_id] = obj.downlink_channel.macroscopic_pathloss_model.cell_assignment(config,obj.id,obj.pos,num_first_UEs);
                if obj.attached_eNodeB.eNodeB_id ~= new_eNodeB_id % 这两个不相等代表要进行切换
                    obj.cell_change.requested = true;
                    obj.cell_change.target_eNodeB = new_eNodeB_id;
                end
            end
        end
        %% 切换函数
        function start_handover(obj,new_eNodeB)
            obj.attached_eNodeB.deattachUser(obj);
            new_eNodeB.attachUser(obj);
        end

        %% 保存每个TTI中UE的SINR
        function save_SINR_in_TTI(obj,config)

            if config.asynchronization && config.isDouble
                obj.ul_dl_SINR_TTI(obj.clock.current_TTI) = obj.ul_dl_wideband_SINR;
                obj.dl_ul_SINR_TTI(obj.clock.current_TTI) = obj.dl_ul_wideband_SINR;
            else
                obj.SINR_TTI(obj.clock.current_TTI) = obj.wideband_SINR;
                obj.ul_SINR_TTI(obj.clock.current_TTI) = obj.ul_wideband_SINR;
            end

            obj.attached_eNodeB_id_TTI(obj.clock.current_TTI) = obj.attached_eNodeB.eNodeB_id;
        end

        %% 下面4个函数用来求SINR
        % 大体思路是先求服务径的CL，再求同频干扰径的CL，然后求邻频干扰径的CL，
        % 在再根据发射功率求接收服务信号功率，干扰信号功率，继而求得SINR

        % 计算上行SINR
        function [IoT interference_level BS_Rx_power BS_Rx_intf_power ul_wideband_loop_SINR] = up_link_quality_model(obj,config,UEs,bs_beam_gain,ue_beam_gain)
            % 输入参数：
            % config：参数配置系统
            % UEs：UE列表
            % bs_beam_gain：BS波束赋型增益
            % ue_beam_gain：UE波束赋型增益
            %
            % 输出参数：
            % IoT：IoT值
            % interference_level：阻塞信号矩阵
            % ul_wideband_loop_SINR：UL上行SINR
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;% 所有基站
            there_are_interferers            = ~isempty(interfering_eNodeBs);

            user_penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; % 得到服务信号的穿损
            user_macroscopic_pathloss = obj.downlink_channel.macroscopic_pathloss + user_penetration_loss; % 得到服务信号的路损（传+单元增益+穿）
            user_shadow_fading_loss = obj.downlink_channel.shadow_fading_pathloss; % 得到服务信号的阴影衰落
            %得到UE的CL
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss-bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)-ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);
            user_CL_linear = 10^(0.1*user_CL);%转换为线性值

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % 计算噪声
            thermal_noise_watts_per_half_RB = obj.BS_thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;

            % 每个UE的最大发射功率
            TX_power = config.UE_tx_power;
            %计算同频干扰及邻频干扰
            if there_are_interferers
                % 求邻频干扰UE（同一个服务基站下的UE干扰）
                interfering_ue_id = [];
                tmp2 = [obj.attached_eNodeB.attached_UEs_vector.id];
                index_in_cell = find(tmp2 == obj.id);% index_in_cell为UE在其所属小区里UE的顺序

                interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;
                for s_ = 1:length(interfering_eNodeBs)
                    tmp = [interfering_eNodeBs(s_).attached_UEs_vector.id];
                    if length(tmp) == length(tmp2)% 本小区与干扰小区UE数相等
                        if isempty(interfering_ue_id)
                            interfering_ue_id = tmp(index_in_cell);
                        else
                            interfering_ue_id = [interfering_ue_id tmp(index_in_cell)];
                        end
                    else% 本小区与干扰小区UE数不相等,实际上不支持每个小区UE数不相同的情形，这里只是提出了一种解决方法

                        if isempty(interfering_ue_id)
                            interfering_ue_id = tmp(NR_randperm(length(tmp),1));
                        else
                            interfering_ue_id = [interfering_ue_id tmp(NR_randperm(length(tmp),1))];
                        end
                    end
                end
                interfering_ue_id = sort(interfering_ue_id); % 为干扰UE排序

                UE_pos_vector = obj.downlink_channel.macroscopic_pathloss_model.UE_pos_vector;
                interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);% 干扰UE的位置

                % 干扰UE到当前UE的服务BS的耦合损耗
                interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);
                interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);
                interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,interfering_ue_id);
                interfering_CL = interfering_macroscopic_pathloss_UE + interfering_penetration_loss + interfering_shadow_fading_loss-bs_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';
                interfering_CL_linear = 10.^(0.1*interfering_CL);

                % 干扰UE到各自服务BS的耦合损耗
                interfering_sector_id = fix((interfering_ue_id-1)/config.UE_per_eNodeB)+1;
                interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(interfering_sector_id,1);
                interfering_enb_id = interfering_enb_id';
                for u_ = 1:length(interfering_ue_id)
                    interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                    interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                    interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),interfering_enb_id(u_));
                    rx_antenna_gain(u_)=ue_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                    beam_antenna_gain(u_)=bs_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                end
                interfering_CL_UE_to_sector = interfering_pathloss_UE_to_sector + interfering_penetrationloss_UE_to_sector + interfering_shadow_fading_UE_to_sector-rx_antenna_gain-beam_antenna_gain;
                interfering_CL_UE_to_sector_linear = 10.^(0.1*interfering_CL_UE_to_sector);

                % 功率控制
                P_min = 10^(-4)/1000;
                R_min = P_min/TX_power;
                CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
                CL_x_ile_linear = 10.^(0.1*CL_x_ile);
                Gama = config.Gama;
                TX_power_PC = TX_power * min(1,max(R_min,(user_CL_linear/CL_x_ile_linear)^Gama));
                interfering_TX_power_PC = TX_power .* min(1,max(R_min,(interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama));
                % 引入功控误差
                TX_power_PC = TX_power_PC * 10.^(0.1*obj.PCe);
                for u_=1:length(interfering_TX_power_PC)
                    interfering_TX_power_PC(u_) = interfering_TX_power_PC(u_) * 10.^(0.1*UEs(interfering_ue_id(u_)).PCe);
                end

                obj.UE_tx_power = TX_power_PC;
                %计算接收功率
                RX_power = TX_power_PC/user_CL_linear;
                %同频干扰接收功率
                interfering_RX_power = interfering_TX_power_PC./interfering_CL_linear;
                %求邻频干扰（双系统下，另一个频率基站的干扰）
                if config.isDouble
                    ACIR_dB = config.ACIR_dB;
                    ACIR_linear = 10^(0.1*ACIR_dB);

                    % 干扰UE
                    % 为应对一个基站同时服务多个UE的问题，这里将一个小区内同时服务的UE分为一组
                    % 如index_in_cell=6，n_UE_served_per_BS=3，则4,5,6为一组
                    if mod(index_in_cell,config.n_UE_served_per_BS) == 0
                        index_in_cell2 = index_in_cell-(config.n_UE_served_per_BS-1):index_in_cell;
                    else
                        yu = mod(index_in_cell,config.n_UE_served_per_BS);% yu=余数，代表该UE在该组内的顺序
                        index_in_cell2 = index_in_cell-(yu-1):index_in_cell+(config.n_UE_served_per_BS-yu);
                    end

                    another_interfering_ue_id = [];
                    another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors;
                    for s_ = 1:length(another_interfering_eNodeBs)
                        tmp = [another_interfering_eNodeBs(s_).attached_UEs_vector.id];
                        if length(tmp) == length(tmp2)% 本小区与干扰小区UE数相等
                            if isempty(another_interfering_ue_id)
                                another_interfering_ue_id = tmp(index_in_cell2);
                            else
                                another_interfering_ue_id = [another_interfering_ue_id tmp(index_in_cell2)];
                            end
                        else% 本小区与干扰小区UE数不相等，实际上不支持每个小区UE数不相同的情形，这里只是提出了一种解决方法
                            if isempty(another_interfering_ue_id)
                                another_interfering_ue_id = tmp(NR_randperm(length(tmp),config.n_UE_served_per_BS));
                            else
                                another_interfering_ue_id = [another_interfering_ue_id tmp(NR_randperm(length(tmp),config.n_UE_served_per_BS))];
                            end
                        end
                    end
                    another_interfering_ue_id = sort(another_interfering_ue_id);
                    another_interfering_ue_pos = UE_pos_vector(another_interfering_ue_id,:);% 干扰UE的位置

                    % 干扰UE到当前UE的服务BS的耦合损耗
                    another_interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,another_interfering_ue_id);
                    another_interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,another_interfering_ue_id);
                    another_interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,another_interfering_ue_id);
                    another_interfering_CL = another_interfering_macroscopic_pathloss_UE + another_interfering_penetration_loss + another_interfering_shadow_fading_loss-bs_beam_gain(another_interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(another_interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';
                    another_interfering_CL_linear = 10.^(0.1*another_interfering_CL);

                    % 干扰UE到各自服务BS的耦合损耗
                    another_interfering_sector_id = fix((another_interfering_ue_id-1)/config.UE_per_eNodeB)+1;
                    another_interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(another_interfering_sector_id,1);
                    another_interfering_enb_id = another_interfering_enb_id';

                    for u_ = 1:length(another_interfering_ue_id)
                        another_interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(another_interfering_ue_pos(u_,:),another_interfering_sector_id(u_));
                        another_interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(another_interfering_ue_pos(u_,:),another_interfering_sector_id(u_));
                        another_interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(another_interfering_ue_pos(u_,:),another_interfering_enb_id(u_));
                        rx_antenna_gain2(u_)=ue_beam_gain(another_interfering_ue_id(u_),UEs(another_interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                        beam_antenna_gain2(u_)=bs_beam_gain(another_interfering_ue_id(u_),UEs(another_interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                    end
                    another_interfering_CL_UE_to_sector = another_interfering_pathloss_UE_to_sector + another_interfering_penetrationloss_UE_to_sector + another_interfering_shadow_fading_UE_to_sector-rx_antenna_gain2-beam_antenna_gain2;
                    another_interfering_CL_UE_to_sector_linear = 10.^(0.1*another_interfering_CL_UE_to_sector);

                    % 功率控制
                    another_interfering_TX_power_PC = TX_power .* min(1,max(R_min,(another_interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama));
                    % 引入功控误差
                    for u_=1:length(another_interfering_TX_power_PC)
                        another_interfering_TX_power_PC(u_) = another_interfering_TX_power_PC(u_) * 10.^(0.1*UEs(another_interfering_ue_id(u_)).PCe);
                    end
                    %邻频干扰接收功率
                    another_interfering_RX_power = another_interfering_TX_power_PC./another_interfering_CL_linear./ACIR_linear;

                    [~,n_i] = size(interfering_RX_power);
                    [~,n_ai] = size(another_interfering_RX_power);
                    if n_ai>n_i
                        % 干扰系统干扰基站数比被干扰系统的多
                        interfering_RX_power(:,(n_i+1:n_ai)) = 0;
                    else
                        another_interfering_RX_power(:,(n_ai+1:n_i)) = 0;
                    end

                    switch config.interference_type
                        case 0
                            total_interfering_power = interfering_RX_power;
                        case 1
                            total_interfering_power = another_interfering_RX_power;
                        case 2
                            total_interfering_power = interfering_RX_power + another_interfering_RX_power;
                    end
                    % 根据ACIR是否使用步进循环进行调整
                    if config.isACIR_loop
                        loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                        loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                        another_interfering_RX_power_1 = another_interfering_RX_power * ACIR_linear; %根据ACIR计算aggressor系统来的干扰强度
                        for i=1:length(loop_ACIR_dB)
                            loop_another_interfering_RX_power(i,:) = another_interfering_RX_power_1 / loop_ACIR_linear(i);
                            loop_total_interfering_power(i,:) = interfering_RX_power + loop_another_interfering_RX_power(i,:);
                            ul_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_total_interfering_power(i,:))+thermal_noise_watts));
                        end
                        %只有同频干扰时的SINR
                        ul_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_RX_power(:))+thermal_noise_watts));
                    else
                        ul_wideband_loop_SINR = 0;
                    end

                    % 阻塞信号
                    another_interfering_RX_power_2 = another_interfering_RX_power * ACIR_linear;
                    interference_level = sum(another_interfering_RX_power_2(:));
                else
                    total_interfering_power = interfering_RX_power;
                    ul_wideband_loop_SINR = 0;
                    interference_level = 0;
                end

                % SINR
                obj.ul_wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))+thermal_noise_watts));
                IoT = 10*log10(sum(total_interfering_power(:))/thermal_noise_watts);% 假如每个小区有10个UE，循环到UE1时，对BS1这个IOT的I就是UE11，21，31的干扰
                % 然后循环UE1到UE10就是BS1取10次快照得到的干扰

                BS_Rx_power = RX_power;
                BS_Rx_intf_power = sum(total_interfering_power(:));
            else
                RX_power = TX_power/user_CL_linear;
                obj.ul_wideband_SINR = 10*log10(RX_power/thermal_noise_watts);
                ul_wideband_loop_SINR = 0;

                BS_Rx_power = RX_power;
                BS_Rx_intf_power = 0;
            end
        end

        function dummy_up_link_quality_model(obj)
            obj.ul_wideband_SINR = NaN;
            obj.UE_tx_power = NaN;
        end

        % 下行SINR
        function [signal_CL interfering_CL another_interfering_CL UE_Rx_power UE_Rx_intf_power wideband_loop_SINR] = down_link_quality_model(obj,config,bs_beam_gain,ue_beam_gain, SYS_config)
            % 输入参数：
            % config：参数配置系统
            % UEs：UE列表
            % bs_beam_gain：BS波束赋型增益
            % ue_beam_gain：UE波束赋型增益
            %
            % 输出参数：
            % signal_CL：服务信号CL
            % interfering_CL：victim系统干扰信号CL
            % another_interfering_CL：aggressor系统干扰信号CL
            % wideband_loop_SINR：UE下行SINR

            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;% 得到victim系统干扰小区
            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors; % 得到aggressor系统干扰小区
            there_are_interferers = ~isempty(interfering_eNodeBs);
            obj.penetration_loss = 0;
%             obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; % 穿损

            d2d = pdist2(obj.attached_site.pos, obj.pos,'Euclidean');
            d3d = sqrt(d2d^2 + (obj.height - obj.attached_eNodeB.tx_height)^2);
            distance = d3d;
            frequency = 4.9e9/1e9;
            pl = 28 + 22*log10(distance) + 20*log10(frequency);
            vertical_angle_grid_el = (180/pi)*(atan2( d2d,obj.height-obj.attached_eNodeB.tx_height));
            vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d,obj.attached_eNodeB.tx_height-obj.height));
            alpha = 90 - vertical_angle_grid_el;
            d3d_cos = d2d/cosd(alpha);
            pl2 = 28 + 22*log10(d3d_cos) + 20*log10(frequency);

            angle_grid = (180/pi)*(atan2((obj.pos(2)-obj.attached_site.pos(2)),(obj.pos(1)-obj.attached_site.pos(1))))-obj.attached_eNodeB.azimuth;
            horizontal_angle_grid = utils.miscUtils.wrapTo359(angle_grid);
            horizontal_angle_grid_s = utils.miscUtils.wrapTo180(horizontal_angle_grid);
            
            angle_grid_ue=180-angle_grid;%UE和BS的水平角互补
            horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(angle_grid_ue);
            horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
            rnd_phi_panel1=180-360*rand(size(obj.attached_eNodeB.tx_height)); %加上-180~180的随机旋转相位
            rnd_phi_panel2=rnd_phi_panel1-180;%两个panel差180度
            phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
            phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
            %加上两侧panel随机旋转角后仍要保证总体方位角在-180到180间
            phi_1 = utils.miscUtils.wrapToAll180(phi_1);
            phi_2 = utils.miscUtils.wrapToAll180(phi_2);


            BS_element_pattern = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
            UE_antenna_gain_1=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
            UE_antenna_gain_2=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
            UE_antenna_gain= max(UE_antenna_gain_1,UE_antenna_gain_2);
            macroscopic_pathloss = pl - BS_element_pattern - UE_antenna_gain;
            macroscopic_pathloss2 = pl2 - BS_element_pattern - UE_antenna_gain;
            user_macroscopic_pathloss        = macroscopic_pathloss ;   % 传输损耗
            user_macroscopic_pathloss_linear = 10^(0.1*user_macroscopic_pathloss);% 转换为线性值
            user_macroscopic_pathloss_linear2 = 10^(0.1*macroscopic_pathloss2);% 转换为线性值
%             user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss;% 阴影衰落
            user_shadow_fading_loss = 0;
            user_shadow_fading_loss_linear   = 10^(0.1*user_shadow_fading_loss); % 转换为线性值

            signal_CL = user_macroscopic_pathloss + user_shadow_fading_loss; % 得到服务信号的CL

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % 计算噪声
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;
            % 发射功率
            TX_power = obj.attached_eNodeB.max_power/config.n_UE_served_per_BS;
            RX_power = TX_power / user_macroscopic_pathloss_linear / user_shadow_fading_loss_linear;
            RX_power2 = TX_power / user_macroscopic_pathloss_linear2;

            if there_are_interferers % 计算干扰信号CL
                parent_sites                            = [interfering_eNodeBs.parent_eNodeB];
                parent_sites_id                         = [parent_sites.id];
                interfering_eNodeB_ids                  = [interfering_eNodeBs.eNodeB_id];
                interfering_RB_grids                    = [interfering_eNodeBs.RB_grid];
                interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
                interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
                interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;% n_干扰enb列，n_RB行
                interfering_power_allocations = sum(interfering_power_allocations);
                for i = 1:length(parent_sites)
                    interfer_sites_pos(i, :) = parent_sites(i).pos;
                    inteNodeB_azimuth(i) = interfering_eNodeBs(i).azimuth;
                end
                d2d_het = pdist2(interfer_sites_pos, obj.pos,'Euclidean');
                d3d_het = sqrt(d2d_het.^2 + (obj.height - obj.attached_eNodeB.tx_height).^2);
                distance_het = d3d_het;
                pl_het = 28 + 22*log10(distance_het) + 20*log10(frequency);
%                 pl_het(1:2) =[150,150];
                vertical_angle_grid_el = (180/pi)*(atan2( d2d_het,obj.height-obj.attached_eNodeB.tx_height));
                vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d_het,obj.attached_eNodeB.tx_height-obj.height));
%                 alpha_het = 90 - vertical_angle_grid_el;
%                 d3d_cos_het = d2d_het./cosd(alpha_het);
                n_i = d2d_het / d2d;

                pl2_het = 28 + 22*log10(d3d_cos) + 22*log10(cosd(alpha)./cosd(alpha./n_i) .* n_i)+ 20*log10(frequency);
                if sum(pl2_het- pl_het)>10
                    pl2_het- pl_het;
                end
                angle_grid = (180/pi)*(atan2((obj.pos(2)-interfer_sites_pos(:,2)),(obj.pos(1)-interfer_sites_pos(:,1))))-inteNodeB_azimuth';
%                 phi = angle_grid;
%                 phi(1) = phi(1)+120;
%                 phi(2) = phi(2)+240;
%                 for i = 3:length(phi)
%                     phi(i) = phi(i) + mod(i,3) * 120;
%                 end
                horizontal_angle_grid = utils.miscUtils.wrapTo359(angle_grid);
                horizontal_angle_grid_s = utils.miscUtils.wrapTo180(horizontal_angle_grid);
                
                angle_grid_ue=180-angle_grid;%UE和BS的水平角互补
                horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(angle_grid_ue);
                horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
                rnd_phi_panel1=180-360*rand(length(parent_sites), 1); %加上-180~180的随机旋转相位
                rnd_phi_panel2=rnd_phi_panel1-180;%两个panel差180度
                phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
                phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
                %加上两侧panel随机旋转角后仍要保证总体方位角在-180到180间
                phi_1 = utils.miscUtils.wrapToAll180(phi_1);
                phi_2 = utils.miscUtils.wrapToAll180(phi_2);


                BS_element_pattern = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
                BS_element_pattern((90-vertical_angle_grid_el)>(SYS_config.tilt+24)) = -inf;
                BS_element_pattern((90-vertical_angle_grid_el)<SYS_config.tilt) = -inf;
                UE_antenna_gain_1=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                UE_antenna_gain_2=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                UE_antenna_gain= max(UE_antenna_gain_1,UE_antenna_gain_2);
                
                macroscopic_pathloss_het = pl_het - BS_element_pattern - UE_antenna_gain;
%                 macroscopic_pathloss_het = macroscopic_pathloss_het(3:3:end);
                macroscopic_pathloss2_het = pl2_het - BS_element_pattern - UE_antenna_gain;
%                 macroscopic_pathloss2_het = real(macroscopic_pathloss2_het);
%                 macroscopic_pathloss2_het = macroscopic_pathloss2_het(3:3:end);

                interfering_macroscopic_pathloss_eNodeB        = macroscopic_pathloss_het ;%干扰路损（求位置时的
%                 interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(parent_sites_id);
                interfering_shadow_fading_loss = 0;
                interfering_CL = interfering_macroscopic_pathloss_eNodeB + interfering_shadow_fading_loss;

                interfering_macroscopic_pathloss_eNodeB_linear = 10.^(0.1*interfering_macroscopic_pathloss_eNodeB);
                interfering_macroscopic_pathloss_eNodeB_linear2 = 10.^(0.1*macroscopic_pathloss2_het);
                interfering_shadow_fading_loss_linear          = 10.^(0.1*interfering_shadow_fading_loss);
                interfering_power = interfering_power_allocations ./ interfering_macroscopic_pathloss_eNodeB_linear';
                interfering_power2 = interfering_power_allocations ./ interfering_macroscopic_pathloss_eNodeB_linear2';
                if config.isDouble %如果是双系统，计算aggressor系统的干扰信号
                    another_parent_sites                            = [another_interfering_eNodeBs.parent_eNodeB];
                    another_parent_sites_id                         = [another_parent_sites.id];
                    another_interfering_eNodeB_ids                  = [another_interfering_eNodeBs.eNodeB_id];
                    another_interfering_RB_grids                    = [another_interfering_eNodeBs.RB_grid];
                    another_interfering_power_allocations_data      = [another_interfering_RB_grids.power_allocation];
                    another_interfering_power_allocations_signaling = [another_interfering_RB_grids.power_allocation_signaling];
                    another_interfering_power_allocations = another_interfering_power_allocations_data + another_interfering_power_allocations_signaling;
                    another_interfering_power = sum(another_interfering_power_allocations);

                    another_interfering_macroscopic_pathloss_eNodeB        = obj.downlink_channel.interfering_macroscopic_pathloss(another_interfering_eNodeB_ids) + obj.downlink_channel.interfering_penetration_loss(another_interfering_eNodeB_ids) - bs_beam_gain(obj.id,another_interfering_eNodeB_ids)'-ue_beam_gain(obj.id,another_interfering_eNodeB_ids)';
                    another_interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(another_parent_sites_id);

                    another_interfering_CL = another_interfering_macroscopic_pathloss_eNodeB + another_interfering_shadow_fading_loss;
                    % 为与同频对比，删掉与UE的服务BS同一个小区的BS
                    [~,index_min] = sort(another_interfering_CL);
                    another_interfering_CL(index_min(1)) = [];

                    another_interfering_macroscopic_pathloss_eNodeB_linear = 10.^(0.1*another_interfering_macroscopic_pathloss_eNodeB);
                    another_interfering_shadow_fading_loss_linear          = 10.^(0.1*another_interfering_shadow_fading_loss);

                    another_interfering_power = another_interfering_power ./ another_interfering_macroscopic_pathloss_eNodeB_linear' ./ another_interfering_shadow_fading_loss_linear';

                    [~,n_i] = size(interfering_power);
                    [~,n_ai] = size(another_interfering_power);
                    if n_ai>n_i
                        % 干扰系统干扰基站数比被干扰系统的多
                        interfering_power(:,(n_i+1:n_ai)) = 0;
                    else
                        another_interfering_power(:,(n_ai+1:n_i)) = 0;
                    end




                    if config.isACIR_loop
                        loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                        loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                        for i=1:length(loop_ACIR_dB)
                            loop_another_interfering_power(i,:) = another_interfering_power / loop_ACIR_linear(i);
                            loop_total_interfering_power(i,:) = interfering_power + loop_another_interfering_power(i,:);
                            wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_total_interfering_power(i,:))+thermal_noise_watts));
                        end
                        wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_power(:))+thermal_noise_watts));
                    else
                        wideband_loop_SINR = 0;
                    end

                    ACIR_dB = config.ACIR_dB;
                    ACIR_linear = 10^(0.1*ACIR_dB);
                    another_interfering_power = another_interfering_power / ACIR_linear;

                    switch config.interference_type
                        case 0
                            total_interfering_power = interfering_power;% 只有本系统干扰
                            total_interfering_power2 = interfering_power2;% 只有本系统干扰
                        case 1
                            total_interfering_power = another_interfering_power;% 只有第二系统干扰
                        case 2
                            total_interfering_power = interfering_power + another_interfering_power; %两个系统的干扰
                    end
                else
                    total_interfering_power = interfering_power;% 只有本系统干扰
                    total_interfering_power2 = interfering_power2;% 只有本系统干扰
                    another_interfering_CL = zeros(size(interfering_power));
                    wideband_loop_SINR = 0;
                end

                % SINR
%                 obj.wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))+thermal_noise_watts));
                % 只有SIR，没有包含噪声N
                obj.wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))));
                obj.wideband_SINR2 = 10*log10(RX_power2/(sum(total_interfering_power(:))+thermal_noise_watts));
%                 obj.wideband_SINR2 = 10*log10(RX_power2/(sum(total_interfering_power2(:))+thermal_noise_watts));
                UE_Rx_power = RX_power;
                UE_Rx_intf_power = sum(total_interfering_power(:));
            else
                wideband_loop_SINR = 0;
                obj.wideband_SINR = 10*log10(RX_power/thermal_noise_watts);

                UE_Rx_power = RX_power;
                UE_Rx_intf_power = 0;
            end

        end

        function dummy_down_link_quality_model(obj)
            obj.wideband_SINR = NaN;
        end

        %% 上行干扰下行（UE干扰UE）只有双系统时才会使用
        function [user_macroscopic_pathloss,ul_down_wideband_loop_SINR] = up_down_model(obj,config,UEs,bs_beam_gain,ue_beam_gain,networkPathlossMap)
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors; %得到干扰eNodeB
            obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; %得到本UE到其服务BS的穿损
            user_macroscopic_pathloss        = obj.downlink_channel.macroscopic_pathloss_900;   % 传+单元增益
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss; %阴影衰落
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss+ obj.penetration_loss -bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)- ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);%计算出CL
            user_CL_linear = 10^(0.1*user_CL);%转换为线性值

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % 噪声
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;
            % 得到基站发射功率
            TX_power = obj.attached_eNodeB.max_power/config.n_UE_served_per_BS;
            % 得到该UE接收到服务BS的功率
            RX_power = TX_power /user_CL_linear;

            %计算出干扰基站的发射功率
            parent_sites                            = [interfering_eNodeBs.parent_eNodeB];
            parent_sites_id                         = [parent_sites.id];
            interfering_eNodeB_ids                  = [interfering_eNodeBs.eNodeB_id];
            interfering_RB_grids                    = [interfering_eNodeBs.RB_grid];
            interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
            interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
            interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;% n_干扰enb列n_RB行
            interfering_power = sum(interfering_power_allocations);

            interfering_macroscopic_pathloss_eNodeB        = obj.downlink_channel.interfering_macroscopic_pathloss(interfering_eNodeB_ids); %传播损耗+单元增益
            interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(parent_sites_id);%阴影衰落
            interfering_CL=interfering_macroscopic_pathloss_eNodeB+interfering_shadow_fading_loss+obj.downlink_channel.interfering_penetration_loss(interfering_eNodeB_ids) - bs_beam_gain(obj.id,interfering_eNodeB_ids)'-ue_beam_gain(obj.id,interfering_eNodeB_ids)';%耦合损耗
            %interfering_PL =    interfering_macroscopic_pathloss_eNodeB-
            interfering_CL_linear = 10.^(0.1*interfering_CL);
            interfering_power = interfering_power ./interfering_CL_linear';%计算接收到的干扰BS功率

            %% 邻系统内采用上行干扰的形式，故需要根据根据干扰UE与被干扰UE计算出两者的单元增益，路径损耗（假设全为NLOS），并且假设无阴影衰落的影响
            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors; %找出邻系统的干扰BS
            another_interfering_eNodeB_ids = [another_interfering_eNodeBs.eNodeB_id]; %找到第二个系统的干扰扇区id号
            ser_group_num = mod(obj.id,config.UE_per_eNodeB); %找出快照组号
            if ser_group_num==0
                ser_group_num=config.UE_per_eNodeB;
            end
            another_interfering_UE_ids = (another_interfering_eNodeB_ids-1).* config.UE_per_eNodeB + ser_group_num;%找出干扰UE的id号
            for u_=1:length(another_interfering_UE_ids)
                %计算该UE被服务时，邻系统小区u_的干扰组
                intf_group =floor((another_interfering_UE_ids(u_)-1)/config.n_UE_served_per_BS)*config.n_UE_served_per_BS+1 :floor((another_interfering_UE_ids(u_)-1)/config.n_UE_served_per_BS)*config.n_UE_served_per_BS+config.n_UE_served_per_BS;
                %计算每个干扰UE相对于该UE的水平角，并转换为(-180,180]
                for j=1:config.n_UE_served_per_BS     %j表示邻系统小区干扰组的第j个UE
                    phi_up_down(u_,j)=atan2(UEs(intf_group(j)).pos(2)-obj.pos(2),UEs(intf_group(j)).pos(1)-obj.pos(1))./pi*180-obj.orientation;
                    phi_up_down(u_,j)=utils.miscUtils.wrapToAll180(phi_up_down(u_,j));
                    %计算每个干扰UE相对于该UE的垂直角
                    theta_up_down(u_,j) = atan2(sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2))^2),UEs(intf_group(j)).height-obj.height)./pi*180;
                    %计算每个干扰UE相对于该UE的单元增益
                    element_gain_interference(u_,j) = obj.antenna.elementPattern(theta_up_down(u_,j),phi_up_down(u_,j));
                    %计算该UE相对于每个干扰UE的垂直角与水平角
                    theta_opposite_observation(u_,j)=atan2(sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2))^2),obj.height-UEs(intf_group(j)).height)./pi*180;
                    phi_opposite_observation(u_,j)=atan2(obj.pos(2)-UEs(intf_group(j)).pos(2),obj.pos(1)-UEs(intf_group(j)).pos(1))./pi*180-UEs(intf_group(j)).orientation;
                    %% 角度调整为(-180,180]
                    phi_opposite_observation(u_,j)=utils.miscUtils.wrapToAll180(phi_opposite_observation(u_,j));

                    %计算该UE相对于每个干扰UE的单元增益
                    element_gain_UE(u_,j) = obj.antenna.elementPattern(theta_opposite_observation(u_,j),phi_opposite_observation(u_,j));
                    %计算二者的3d距离
                    d_3d(u_,j) = sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2)).^2+(obj.height-UEs(intf_group(j)).height).^2);
                    %计算pathloss
                    if obj.id <= networkPathlossMap.num_first_UEs
                        frequency=config.frequency2;
                    else
                        frequency=config.frequency;
                    end
                    switch config.scene_type %均采用NLOS路损模型
                        case {'InH','InH2','InH3'}%如果是indoor场景则采用室内路损模型
                            shadow_fading_std = 8.03;%采用室内NLOS时的阴影衰落
                            %                             shadow_fading_std = 3.0;%采用室内LOS时的阴影衰落
                            shadow_fading_dB(u_,j)=normrnd(0,shadow_fading_std);%产生对数正态分布的阴影衰落
                            pl(u_,j)=17.3+38.3*log10(d_3d(u_,j))+24.9*log10(frequency/1e9)+shadow_fading_dB(u_,j);%NLOS
                            %                             pl(u_,j)=32.4+17.3*log10(d_3d(u_))+20*log10(frequency/1e9);%LOS
                        otherwise %XIA模型中采用 天线低于平均房顶水平的模型,包括阴影，传播损耗和绕射损耗。“A Simplified Analytical Model for Predicting Path Loss in Urban and Suburban Environments”
                            h_roof = 12;%屋顶高度12m
                            w=30;%街道宽度30m
                            d=80;%两建筑中心距离
                            R=d_3d(u_,j); %发射机与接收机间距
                            delta_hb= obj.height - h_roof;%接收机相对于房顶的高度（比房顶低则为负）
                            delta_hm = h_roof - UEs(intf_group(j)).height;%发射机相对于房顶的高度差
                            phi=-atan2(delta_hb,d);%到第一排建筑物边缘的入射角
                            lamda = 3e8/config.frequency; %波长
                            x=w/2;%干扰UE相对于衍射边的水平距离
                            theta=atan2(delta_hm,x);
                            r=sqrt(delta_hm^2+x^2);%
                            pl(u_,j)=-10*log10((lamda/(2*sqrt(2)*pi*R))^2)-10*log10(lamda/(2*pi^2*r)*((1/theta)-1/(2*pi+theta))^2)-10*log((d/(2*pi*(R-d)))^2*lamda/sqrt(delta_hb^2+d^2)*(1/phi-1/(2*pi+phi))^2);
                    end

                    %计算穿透损耗
                    switch config.scene_type
                        case {'UMA','UMI','RMa','UMa_to_UMi'}
                            penetration(u_,j)=obj.downlink_channel.macroscopic_penetration_loss+UEs(intf_group(j)).downlink_channel.macroscopic_penetration_loss;
                        case {'InH','InH2'}
                            penetration(u_,j)=0;
                        case {'UMi_to_InH','UMa_to_InH'} %室内相对于在同一室内的用户无穿透损耗，但相对于Urban场景下的用户一定会穿一堵墙
                            L_glass = 2+0.2*config.frequency/1e9;
                            L_concrete = 5+4*config.frequency/1e9;
                            N_low = 0+4.4*randn(1);
                            low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                            L_IRRglass = 23+0.3*config.frequency/1e9;
                            N_high = 0+6.5*randn(1);
                            high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                            if rand(1)>0.5
                                penetration(u_,j)=low;
                            else
                                penetration(u_,j)=high;
                            end
                        case {'InH3'} %室内相对于在同一室内的用户无穿透损耗，但相对于Urban场景下的用户一定会穿一堵墙
                            L_glass = 2+0.2*config.frequency/1e9;
                            L_concrete = 5+4*config.frequency/1e9;
                            N_low = 0+4.4*randn(1);
                            low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                            L_IRRglass = 23+0.3*config.frequency/1e9;
                            N_high = 0+6.5*randn(1);
                            high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                            if rand(1)>0.5
                                penetration(u_,j)=2*low;
                            else
                                penetration(u_,j)=2*high;
                            end

                    end
                    %计算耦合损耗

                    another_interfering_CL(u_,j) = pl(u_,j) + penetration(u_,j)-element_gain_interference(u_,j) - element_gain_UE(u_,j) - ue_beam_gain(obj.id,UEs(intf_group(j)).attached_eNodeB.eNodeB_id) - bs_beam_gain(obj.id,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);
                    another_interfering_CL_linear(u_,j) = 10.^(0.1*another_interfering_CL(u_,j)); %转化为线性值
                    %计算干扰UE的功控后发射功率（通过UE对服务基站的CL来确定）
                    another_interfering_pathloss_UE_to_sector(u_,j) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);%干扰UE到其服务小区的传播损耗
                    another_interfering_penetrationloss_UE_to_sector(u_,j) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);%干扰UE到其服务小区的穿透损耗
                    another_interfering_shadow_fading_UE_to_sector(u_,j) = obj.downlink_channel.shadow_fading_model.get_pathloss(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.parent_eNodeB.id);%干扰UE到其服务小区的阴影衰落
                    tx_antenna_gain(u_,j)=ue_beam_gain(intf_group(j),UEs(intf_group(j)).attached_eNodeB.eNodeB_id); %干扰UE到其服务小区的UE侧波束赋型增益
                    rx_antenna_gain(u_,j)=bs_beam_gain(intf_group(j),UEs(intf_group(j)).attached_eNodeB.eNodeB_id); %干扰UE到其服务小区的BS侧波束赋型增益


                    %干扰UE到其服务小区的耦合损耗（用于计算功控）
                    another_interfering_CL_UE_to_sector(u_,j) = another_interfering_pathloss_UE_to_sector(u_,j) + another_interfering_penetrationloss_UE_to_sector(u_,j) + another_interfering_shadow_fading_UE_to_sector(u_,j)-rx_antenna_gain(u_,j)-tx_antenna_gain(u_,j);
                    another_interfering_CL_UE_to_sector_linear(u_,j) = 10.^(0.1*another_interfering_CL_UE_to_sector(u_,j));
                    another_UE_up(u_,j) = element_gain_UE(u_,j)+rx_antenna_gain(u_,j);

                    another_BS_up(u_,j) =  element_gain_interference(u_,j);%+tx_antenna_gain(u_,j);



                    % 每个UE的最大发射功率
                    TX_power_ue = config.UE_tx_power;
                    %功率控制
                    P_min = 10^(-4)/1000;
                    R_min = P_min/TX_power_ue;
                    CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
                    CL_x_ile_linear = 10.^(0.1*CL_x_ile);
                    Gama = config.Gama;
                    %得到干扰UE功控后的发射
                    another_interfering_TX_power(u_,j) = TX_power_ue * min(1,max(R_min,(another_interfering_CL_UE_to_sector_linear(u_,j)/CL_x_ile_linear).^Gama));
                end
            end

            obj.ul_ad_TX_power = another_interfering_TX_power;
            obj.ul_dl_UE_gain = another_UE_up;
            obj.ul_dl_BS_gain = another_BS_up;
            obj.ue_ue_pathloss = pl+penetration;
            if config.isACIR_loop %循环记录ACIR
                loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                for i=1:length(loop_ACIR_dB)
                    loop_another_interfering_power(i,:) = sum(another_interfering_TX_power./another_interfering_CL_linear./ loop_ACIR_linear(i),2)';%得到接收到邻系统干扰UE组的干扰功率之和
                    ul_down_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_another_interfering_power(i,:))+sum(interfering_power(:))+thermal_noise_watts));%计算上行干扰下行的SINR
                end
                ul_down_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_power(:))+thermal_noise_watts));%计算出ACIR无穷大时的SINR
            else
                ul_down_wideband_loop_SINR = 0;
            end
            ACIR_dB = config.ACIR_dB;%用于校准
            ACIR_linear = 10^(0.1*ACIR_dB);
            another_interfering_RX_power = sum(another_interfering_TX_power./another_interfering_CL_linear./ACIR_linear,2)'; %得到接收到邻系统干扰UE组的干扰功率之和
            % Interfering_power_asy = interfering_power(:)+another_interfering_RX_power(:);
            total_interfering_power=sum(interfering_power(:))+sum(another_interfering_RX_power(:));
            obj.ul_dl_wideband_SINR=10*log10(RX_power/(total_interfering_power+thermal_noise_watts));

        end

        %%  下行干扰上行（BS干扰BS）只有双系统时才会使用
        function [user_macroscopic_pathloss,down_ul_wideband_loop_SINR] = down_up_model(obj,config,UEs,bs_beam_gain,ue_beam_gain,networkPathlossMap)
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors; %得到干扰eNodeB
            obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; %得到该UE到其服务BS的穿损
            user_macroscopic_pathloss        = obj.downlink_channel.macroscopic_pathloss;%_900;   %传+单元增益
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss;%得到该UE到其服务BS的阴影衰落
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss + obj.penetration_loss -bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)- ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);%计算出CL
            user_CL_linear = 10^(0.1*user_CL);%转换为线性值

            %BS_antenna_gain = obj.downlink_channel.macroscopic_pathloss_model.BS_antenna_gain;
            % BS_antenna_gain = networkPathlossMap.BS_antenna_gain;
            % UE_antenna_gain = networkPathlossMap.UE_antenna_gain;

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % 噪声
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;

            interfering_ue_id = [];
            tmp2 = [obj.attached_eNodeB.attached_UEs_vector.id]; %得到该UE服务BS的服务用户集
            index_in_cell = find(tmp2 == obj.id);% index_in_cell为UE在其所属小区里UE的索引

            for s_ = 1:length(interfering_eNodeBs) %得到干扰基站正在服务的用户集
                tmp = [interfering_eNodeBs(s_).attached_UEs_vector.id];%得到干扰基站的服务用户集
                if length(tmp) == length(tmp2)
                    if isempty(interfering_ue_id) %得到干扰UE
                        interfering_ue_id = tmp(index_in_cell);
                    else
                        interfering_ue_id = [interfering_ue_id tmp(index_in_cell)];
                    end
                else

                    if isempty(interfering_ue_id)
                        interfering_ue_id = tmp(NR_randperm(length(tmp),1));
                    else
                        interfering_ue_id = [interfering_ue_id tmp(NR_randperm(length(tmp),1))];
                    end
                end
            end
            interfering_ue_id = sort(interfering_ue_id); % 为干扰UE排序

            UE_pos_vector = obj.downlink_channel.macroscopic_pathloss_model.UE_pos_vector;
            interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);% 干扰UE的位置
            % 干扰UE到当前UE的服务BS的耦合损耗（本系统是上行干扰）
            interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);%传+单元
            interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);%穿
            interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,interfering_ue_id);%阴衰
            interfering_CL = interfering_macroscopic_pathloss_UE + interfering_penetration_loss + interfering_shadow_fading_loss-bs_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';%耦合损耗
            interfering_CL_linear = 10.^(0.1*interfering_CL);

            % 干扰UE到各自服务BS的耦合损耗（用于计算干扰UE的功控）
            interfering_sector_id = fix((interfering_ue_id-1)/config.UE_per_eNodeB)+1; %干扰扇区号
            interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(interfering_sector_id,1); %干扰eNodeB号
            interfering_enb_id = interfering_enb_id';

            for u_ = 1:length(interfering_ue_id)
                interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));%传+单元增益

                interfering_UE_gain(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_ue_antenna_gain(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                interfering_BS_gain(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_bs_antenna_gain(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));%穿透损耗
                interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),interfering_enb_id(u_));%阴影衰落
                rx_antenna_gain(u_)=ue_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);%UE端波束赋型增益
                beam_antenna_gain(u_)=bs_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);%BS端波束赋型增益

                interfering_gain_UE(u_) = interfering_UE_gain(u_); %+ rx_antenna_gain(u_);
                interfering_gain_BS(u_)  = interfering_BS_gain(u_) + beam_antenna_gain(u_);
            end
            interfering_CL_UE_to_sector = interfering_pathloss_UE_to_sector + interfering_penetrationloss_UE_to_sector + interfering_shadow_fading_UE_to_sector-rx_antenna_gain-beam_antenna_gain;%干扰UE到其服务BS的耦合损耗
            interfering_CL_UE_to_sector_linear = 10.^(0.1*interfering_CL_UE_to_sector);

            obj.interfering_BS = interfering_gain_BS;
            obj.interfering_UE = interfering_gain_UE;

            % 每个UE最大发射功率
            TX_power_ue = config.UE_tx_power;
            % 功率控制
            P_min = 10^(-4)/1000;
            R_min = P_min/TX_power_ue;
            CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
            CL_x_ile_linear = 10.^(0.1*CL_x_ile);
            Gama = config.Gama;
            TX_power_PC = TX_power_ue * min(1,max(R_min,(user_CL_linear/CL_x_ile_linear)^Gama));%该UE的发射功率

            % 为该UE引入功控误差
            TX_power_PC = TX_power_PC * 10.^(0.1*obj.PCe);
            interfering_TX_power_PC = TX_power_ue .* min(1,max(R_min,(interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama)); %干扰UE的发射功率
            for u_=1:length(interfering_TX_power_PC)
                interfering_TX_power_PC(u_) = interfering_TX_power_PC(u_) * 10.^(0.1*UEs(interfering_ue_id(u_)).PCe);%为干扰UE引入功控误差
            end

            obj.UE_tx_power = TX_power_PC;
            RX_power = TX_power_PC/user_CL_linear; %求出该UE服务BS的接收功率
            interfering_RX_power = interfering_TX_power_PC./interfering_CL_linear; %求出接收到同一系统下其他干扰UE的功率

            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors;%得到邻系统干扰扇区
            another_interfering_eNodeB_ids = [another_interfering_eNodeBs.eNodeB_id]; %找到第二个系统的干扰扇区id号
            for enb=1:length(another_interfering_eNodeB_ids)
                %计算每个干扰BS观测到服务BS的水平角，并转换为(-180,180]
                phi_down_up(enb)=atan2(another_interfering_eNodeBs(enb).parent_eNodeB.pos(2)-obj.attached_eNodeB.parent_eNodeB.pos(2),another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)-obj.attached_eNodeB.parent_eNodeB.pos(1))./pi*180-obj.attached_eNodeB.azimuth;
                phi_down_up(enb)=utils.miscUtils.wrapToAll180(phi_down_up(enb));
                %计算每个干扰UE到该UE的垂直角
                switch config.scene_type
                    case {'RMa','UMA','UMI','InH','InH2','InH3'} %同构下双系统站址高度相同
                        h1(enb)=config.site_height; %h1为victim的站高
                        h2(enb)=config.site_height; % h2为aggressor系统的站高
                        another_tx_power(enb)=config.eNodeB_tx_power/config.n_UE_served_per_BS;%发射功率与带宽有关，继而与同时服务的用户数有关
                    otherwise %异构下站址高度不同
                        if another_interfering_eNodeB_ids(enb)> networkPathlossMap.num_first_sectors %干扰基站是第二个系统的(高度差为第二个系统站高于第一个系统的站高差值)
                            h1(enb)=config.site_height;
                            h2(enb)=config.site_height2;
                            another_tx_power(enb)=config.eNodeB_tx_power2/config.n_UE_served_per_BS;
                        else %干扰基站是第一个系统的(高度差为第一个系统站高于第二个系统的站高差值)
                            h1(enb)=config.site_height2;
                            h2(enb)=config.site_height;
                            another_tx_power(enb)=config.eNodeB_tx_power/config.n_UE_served_per_BS;
                        end
                end
                %干扰BS侧观测到服务BS的垂直角
                theta_down_up(enb) = atan2(sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2))^2), h2(enb)- h1(enb))./pi*180;
                %计算每个干扰BS相对该UE所归属的BS的单元增益
                element_gain_interference(enb) = another_interfering_eNodeBs(enb).antenna.elementPattern(theta_down_up(enb),phi_down_up(enb));
                %计算该UE所归属的BS相对于每个干扰BS的单元增益
                theta_opposite_observation(enb)=atan2(sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2))^2), h1(enb)- h2(enb))./pi*180;
                phi_opposite_observation(enb)=atan2(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2),obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1))./pi*180-another_interfering_eNodeBs(enb).azimuth;

                %% 角度调整为(-180,180]
                phi_opposite_observation(enb)=utils.miscUtils.wrapToAll180(phi_opposite_observation(enb));
                %计算victim系统基站相对于干扰BS的单元增益
                element_gain_victim(enb) = obj.attached_eNodeB.antenna.elementPattern(theta_opposite_observation(enb),phi_opposite_observation(enb));

                %计算二者的3d距离
                d_3d(enb) = sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2)).^2+(h2(enb)-h1(enb)).^2);
                if obj.id <= networkPathlossMap.num_first_UEs
                    frequency=config.frequency2;
                else
                    frequency=config.frequency;
                end
                if strcmp(config.macroscopic_pathloss_model,'TS38900') && strcmp(config.macroscopic_pathloss_model2,'TS38900') %如果两个系统都是NR，则用NR的计算公式
                    switch config.scene_type %均采用NLOS路损模型，异构时采用大场景的传播损耗公式
                        case {'UMA','UMa_to_InH','UMa_to_UMi'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'urban_macro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMA_std = 6;%UMA NLOS 的阴影衰落标准差
                            else
                                UMA_std = 4;%UMA LOS 的阴影衰落标准差
                            end
                            sf_UMA(enb) = normrnd(0,UMA_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_UMA(enb);
                            %                               pl(enb)=13.54+39.08*log10(d_3d(enb))+20*log10(frequency/1e9)-0.6*(obj.height-1.5)+sf_UMA(enb);
                        case {'UMI','UMi_to_InH'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'urban_micro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMI_std = 7.82;%UMI NLOS 的阴影衰落标准差
                            else
                                UMI_std = 4;%UMI LOS 的阴影衰落标准差
                            end
                            sf_UMI(enb) = normrnd(0,UMI_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_UMI(enb);
                            %                               pl(enb)=35.3*log10(d_3d(enb))+22.4+21.3*log10(frequency/1e9)-0.3*(obj.height-1.5)+sf_UMI(enb);
                        case {'InH','InH2','InH3'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'indoor',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                InH_std = 8.03;%UMA NLOS 的阴影衰落标准差
                            else
                                InH_std = 3;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_InH(enb) = normrnd(0,InH_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_InH(enb);
                            %                               pl(enb)=17.3+38.3*log10(d_3d(enb))+24.9*log10(frequency/1e9)+sf_InH(enb);%NLOS
                            %                           pl(enb)=32.4+17.3*log10(d_3d(enb))+20*log10(frequency/1e9);%LOS
                        case 'RMa'
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'rural_macro',config);
                            temp_model.h_ut = 35;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                RMA_std = 8;%UMA NLOS 的阴影衰落标准差
                            else
                                RMA_std = 4;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_RMA(enb) = normrnd(0,RMA_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_RMA(enb);
                            %                               pl(enb)=161.04-7.1*log10(20)+7.5*log10(5)-(24.37-3.7*(5/35).^2)*log10(35)+(43.42-3.1*log10(35))*log10(d_3d(enb)-3)+20*log10(frequency/1e9)-(3.2*(log10(11.75*obj.height)).^2-4.97)+sf_RMa(enb);
                    end
                elseif strcmp(config.macroscopic_pathloss_model,'TS38901') && strcmp(config.macroscopic_pathloss_model2,'TS38901')
                    switch config.scene_type %均采用NLOS路损模型，异构时采用大场景的传播损耗公式
                        case {'UMA','UMa_to_InH','UMa_to_UMi'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'urban_macro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMA_std = 6;%UMA NLOS 的阴影衰落标准差
                            else
                                UMA_std = 4;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_UMA(enb) = normrnd(0,UMA_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_UMA(enb);
                            %                               pl(enb)=13.54+39.08*log10(d_3d(enb))+20*log10(frequency/1e9)-0.6*(obj.height-1.5)+sf_UMA(enb);
                        case {'UMI','UMi_to_InH'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'urban_micro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMI_std = 7.82;%UMA NLOS 的阴影衰落标准差
                            else
                                UMI_std = 4;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_UMI(enb) = normrnd(0,UMI_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_UMI(enb);
                            %                               pl(enb)=35.3*log10(d_3d(enb))+22.4+21.3*log10(frequency/1e9)-0.3*(obj.height-1.5)+sf_UMI(enb);
                        case {'InH','InH2','InH3'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'indoor',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                InH_std = 8.03;%UMA NLOS 的阴影衰落标准差
                            else
                                InH_std = 3;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_InH(enb) = normrnd(0,InH_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_InH(enb);
                            %                               pl(enb)=17.3+38.3*log10(d_3d(enb))+24.9*log10(frequency/1e9)+sf_InH(enb);%NLOS
                            %                           pl(enb)=32.4+17.3*log10(d_3d(enb))+20*log10(frequency/1e9);%LOS
                        case 'RMa'
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'rural_macro',config);
                            temp_model.h_ut = 35;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                RMA_std = 8;%UMA NLOS 的阴影衰落标准差
                            else
                                RMA_std = 4;%UMA NLOS 的阴影衰落标准差
                            end
                            sf_RMA(enb) = normrnd(0,RMA_std); %对数正态分布
                            pl(enb) = pl(enb) + sf_RMA(enb);
                            %                               pl(enb)=161.04-7.1*log10(20)+7.5*log10(5)-(24.37-3.7*(5/35).^2)*log10(35)+(43.42-3.1*log10(35))*log10(d_3d(enb)-3)+20*log10(frequency/1e9)-(3.2*(log10(11.75*obj.height)).^2-4.97)+sf_RMa(enb);
                    end
                else %UMA或RMA场景为LTE
                    LTE_std=10;%LTE下阴影衰落标准差为10dB
                    LTE_sf(enb) = normrnd(0,LTE_std); %对数正态分布
                    switch config.scene_type
                        case 'RMa'%如果是RMa场景
                            pl(enb) = 69.55+26.16*log10(config.frequency/1e6)-13.82*log10(config.site_height)+(44.9-6.55*log10(config.site_height))*log10(d_3d(enb)/1e3)-4.78*(log10(config.frequency/1e6))^2+18.33*log10(config.frequency/1e6)-40.94+LTE_sf(enb);
                        otherwise %如果是UMa场景或者包含UMA的异构场景
                            pl(enb) = 40*(1-4e-3*abs(h2(enb)-h1(enb)))*log10(d_3d(enb)/1e3)-18*log10(abs(h2(enb)-h1(enb)))+21*log10(config.frequency/1e6)+80+LTE_sf(enb);
                    end

                end
                switch config.scene_type
                    case {'UMA','UMI','RMa','InH','UMa_to_UMi','InH2'}
                        penetration(enb)=0;%基站干扰基站无穿透损耗
                    case {'UMi_to_InH','UMa_to_InH'} %室外对室内的异构场景考虑穿透一堵墙
                        L_glass = 2+0.2*config.frequency/1e9;
                        L_concrete = 5+4*config.frequency/1e9;
                        N_low = 0+4.4*randn(1);
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                        L_IRRglass = 23+0.3*config.frequency/1e9;
                        N_high = 0+6.5*randn(1);
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                        if rand(1)>0.5 % 50%低穿损，50%高穿损
                            penetration(enb)=low;
                        else
                            penetration(enb)=high;
                        end
                    case {'InH3'} %室外对室内的异构场景考虑穿透一堵墙
                        L_glass = 2+0.2*config.frequency/1e9;
                        L_concrete = 5+4*config.frequency/1e9;
                        N_low = 0+4.4*randn(1);
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                        L_IRRglass = 23+0.3*config.frequency/1e9;
                        N_high = 0+6.5*randn(1);
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                        if rand(1)>0.5 % 50%低穿损，50%高穿损
                            penetration(enb)=low*2;
                        else
                            penetration(enb)=high*2;
                        end
                end
                %计算邻系统干扰BS到本UE服务BS的耦合损耗
                another_interfering_CL(enb) = pl(enb) + penetration(enb)-element_gain_interference(enb) - element_gain_victim(enb) - ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) - bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_sent(enb) = pl(enb) + penetration(enb)-element_gain_victim(enb)  - ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) ;
                another_UE(enb) = -element_gain_interference(enb)-ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) ;
                another_BS(enb) =  element_gain_victim(enb)+bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_inter(enb) = pl(enb) + penetration(enb)-element_gain_interference(enb)  - bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_linear(enb) = 10.^(0.1*another_interfering_CL(enb)); %转化为线性值
            end
            %               obj.bs_bs_pathloss = pl + penetration;%路损
            obj.dl_ul_UE_gain = another_UE;
            obj.dl_ul_BS_gain = another_BS;
            %               obj.bs_bs_pathloss = another_interfering_CL_sent;%路损+发射天线增益
            obj.bs_bs_pathloss = another_interfering_CL_inter ;%路损+接收天线增益

            if config.isACIR_loop %是否循环统计不同ACIR下的结果
                loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                for i=1:length(loop_ACIR_dB)
                    loop_another_interfering_power(i,:) = another_tx_power./another_interfering_CL_linear./ loop_ACIR_linear(i);%计算收到的邻系统干扰BS泄露的功率
                    %干扰基站对共址的邻系统基站的干扰按路损公式算无限大，因而假设共址的两个系统基站完全隔离
                    for j=1:length(loop_another_interfering_power(i,:))
                        if loop_another_interfering_power(i,j)==Inf || isnan(loop_another_interfering_power(i,j))
                            loop_another_interfering_power(i,j)=0; %邻系统共址BS对本系统的BS的干扰设置为0
                        end
                    end
                    down_ul_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_another_interfering_power(i,:))+sum(interfering_RX_power(:))+thermal_noise_watts));%计算下行干扰上行的SINR
                end
                down_ul_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_RX_power(:))+thermal_noise_watts));
            else
                down_ul_wideband_loop_SINR = 0;
            end
            ACIR_dB = config.ACIR_dB; %校准使用的ACIR值
            ACIR_linear = 10^(0.1*ACIR_dB);
            another_interfering_RX_power=another_tx_power./another_interfering_CL_linear./ACIR_linear; %得到接收到的干扰UE功率;
            %干扰基站对共址的邻系统基站的干扰按路损公式算无限大，因而假设共址的两个系统基站完全隔离
            total_interfering_power=sum(interfering_RX_power(:))+sum(another_interfering_RX_power(find(another_interfering_RX_power~=inf & ~isnan(another_interfering_RX_power))));
            obj.dl_ul_wideband_SINR=10*log10(RX_power/(total_interfering_power+thermal_noise_watts));  %将校准时使用的ACIR下的SINR记录下来
        end

        function dummy_up_down_model(obj)
            obj.ul_dl_wideband_SINR = NaN;
        end

        function dummy_down_up_model(obj)
            obj.dl_ul_wideband_SINR = NaN;
        end

        function struct_out = basic_information_in_struct(obj)
            struct_out.id                    = obj.id;
            struct_out.pos                   = obj.pos;
        end

        function distance_to_site = distance_to_attached_site(obj)
            % Return the distance to the attached site (function added for convenience)
            distance_to_site = sqrt(sum((obj.attached_site.pos -obj.pos).^2));
        end

        % Clear all non-basic info and leaves just basic information describing the UE
        function clear_non_basic_info(obj)
            obj.attached_site             = [];
            obj.attached_sector_idx       = [];
            obj.attached_eNodeB           = [];
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.clock                     = [];
            obj.wideband_SINR             = [];
            obj.ul_wideband_SINR          = [];
            obj.ul_dl_wideband_SINR       = [];
            obj.dl_ul_wideband_SINR       = [];
            obj.SINR_TTI                  = [];
            obj.ul_SINR_TTI               = [];
            obj.ul_dl_SINR_TTI            = [];
            obj.dl_ul_SINR_TTI            = [];
            obj.attached_eNodeB_id_TTI    = [];
            obj.deactivate_UE             = [];
            obj.receiver_noise_figure     = [];
            obj.thermal_noise_W_RB        = [];
            obj.penetration_loss          = [];
            obj.rx_antenna_gain           = [];
            obj.beam_antenna_gain         = [];
        end

    end
end
