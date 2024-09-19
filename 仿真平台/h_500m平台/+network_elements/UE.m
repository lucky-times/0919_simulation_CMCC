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
            rng(0);
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;% 得到victim系统干扰小区
            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors; % 得到aggressor系统干扰小区
            there_are_interferers = ~isempty(interfering_eNodeBs);
            obj.penetration_loss = 0;
%             obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; % 穿损

            d2d = pdist2(obj.attached_site.pos, obj.pos,'Euclidean');
            d3d = sqrt(d2d^2 + (obj.height - obj.attached_eNodeB.tx_height)^2);
            distance = d3d;
            frequency = config.frequency;
            pl = NR_FSPL(distance, frequency);
%             pl = 32.45 + 22*log10(distance) + 20*log10(frequency);
            vertical_angle_grid_el = (180/pi)*(atan2( d2d,obj.height-obj.attached_eNodeB.tx_height));
            vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d,obj.attached_eNodeB.tx_height-obj.height));
            alpha = 90 - vertical_angle_grid_el;
            d3d_cos = d2d/cosd(alpha);
            pl2 = NR_FSPL(d3d_cos, frequency);
%             pl2 = 32.45 + 22*log10(d3d_cos) + 20*log10(frequency);
            ang = (180/pi)*(atan2((obj.pos(2)-obj.attached_site.pos(2)),(obj.pos(1)-obj.attached_site.pos(1))));
            angle_grid = mod(ang, 360)-obj.attached_eNodeB.azimuth;
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
%                 interfering_eNodeB_height               = [interfering_eNodeBs.heighjt];
                interfering_RB_grids                    = [interfering_eNodeBs.RB_grid];
                interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
                interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
                interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;% n_干扰enb列，n_RB行
                interfering_power_allocations = sum(interfering_power_allocations);
                for i = 1:length(parent_sites)
                    interfer_sites_pos(i, :) = parent_sites(i).pos;
                    inteNodeB_azimuth(i) = interfering_eNodeBs(i).azimuth;
                    interfering_eNodeB_height(i,1) = interfering_eNodeBs(i).tx_height;
                end
                d2d_het = pdist2(interfer_sites_pos, obj.pos,'Euclidean');
                d3d_het = sqrt(d2d_het.^2 + (obj.height - obj.attached_eNodeB.tx_height).^2);
                distance_het = d3d_het;
                %由三维距离计算
                pl_het =NR_FSPL(distance_het,frequency);
                
%                 pl_het(1:2) =[150,150];
                vertical_angle_grid_el = (180/pi)*(atan2( d2d_het,obj.height-interfering_eNodeB_height));
                vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d_het,interfering_eNodeB_height-obj.height));
%                 alpha_het = 90 - vertical_angle_grid_el;
%                 d3d_cos_het = d2d_het./cosd(alpha_het);
                n_i = d2d_het / d2d;
                %按照论文所提公式计算
                pl2_het = NR_FSPL(d3d_cos,frequency) + 22*log10(cosd(alpha)./cosd(alpha./n_i) .* n_i);
                if sum(pl2_het- pl_het)>0
                    tmp = pl2_het - pl_het;
                end
                for i = 1:length(pl2_het)
                    if isreal(pl2_het(i))
                        continue;
                    end
                    pl2_het(i) = inf;
                end
                ang = (180/pi)*(atan2((obj.pos(2)-interfer_sites_pos(:,2)),(obj.pos(1)-interfer_sites_pos(:,1))));
                angle_grid = mod(ang, 360)-inteNodeB_azimuth';
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
                beta = 10^-0.3;%干扰小区天线的缩减因子
                reduction_factor = 10*log10(beta);
                BS_element_pattern = zeros(size(vertical_angle_grid_el));
                for i = 1 : length(vertical_angle_grid_el)
                    BS_element_pattern(i) = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el(i), horizontal_angle_grid_s(i), SYS_config.tilt) + reduction_factor ;%均值为-10，方差为10的随机数，正态分布   + sqrt(10)*randn()-10
                end
%                 BS_element_pattern = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
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

                % SINR,加扰50%，干扰功率将为原先的一半
                obj.wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))/2+thermal_noise_watts));
                obj.wideband_SINR2 = 10*log10(RX_power2/(sum(total_interfering_power2(:))/2+thermal_noise_watts));
                UE_Rx_power = RX_power;
                UE_Rx_intf_power = sum(total_interfering_power(:));
            else
                wideband_loop_SINR = 0;
                obj.wideband_SINR = 10*log10(RX_power/thermal_noise_watts);

                UE_Rx_power = RX_power;
                UE_Rx_intf_power = 0;
            end

        end

       

    end
end
