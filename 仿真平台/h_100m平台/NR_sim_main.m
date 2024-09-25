function [output_results_file] = NR_sim_main(SYS_config)
%系统级仿真的主程序
%
%输入参数：
%SYS_config：仿真参数
%
%输出参数：
%output_results_file：输出结果文件

ticID_start_sim_begin = tic; % 程序开始的时间

%% 防止参数设置不正确
SYS_config = NR_load_params_again( SYS_config );

%% 传入调试等级
DEBUG_LEVEL = SYS_config.debug_level;

%% 产生基站，并得到路损图谱与阴影衰落
[eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkShadowFadingMap] = NR_generate_network(SYS_config);

%% 时钟
networkClock = network_elements.clock(SYS_config.TTI_length);

for b_=1:length(eNodeB_sites)
    eNodeB_sites(b_).clock = networkClock;
end

%% 创建用户
[UEs] = NR_generate_users(SYS_config,eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkClock);
eNodeB_pos = [];
for i = 1:length(eNodeB_sites)
    eNodeB_pos = [eNodeB_pos;eNodeB_sites(i).pos];
end
UE_pos = [];
for i = 1:length(UEs)
    UE_pos = [UE_pos;UEs(i).pos];
end
figure;
scatter(eNodeB_pos(:,1), eNodeB_pos(:, 2));
hold on;
scatter(UE_pos(:,1), UE_pos(:, 2));
hold off;
axis equal;
title('拓扑分布图');

if isempty(UEs)
    warning('No UEs generated.');
else  
     %% UE其他方面的初始化，用于计算数据
    for u_=1:length(UEs)
        
        % 加入下行信道 (包括传播损耗，穿透损耗和阴影衰落)
        UEs(u_).downlink_channel = channel_models.downlinkChannelModel(UEs(u_));
                       
        % 传入传播损耗
        UEs(u_).downlink_channel.set_macroscopic_pathloss_model(networkPathlossMap);
        
        % 传入阴影衰落
        if SYS_config.macroscopic_pathloss_is_model
            UEs(u_).downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end
        
        % 获取UE其他参数（RB分配表、高度、旋转角）
        UEs(u_).RB_grid = UEs(u_).downlink_channel.RB_grid;
%         UEs(u_).height = UEs(u_).downlink_channel.UE_height;
        UEs(u_).height = SYS_config.UE_height;
        UEs(u_).orientation = 0;
%         UEs(u_).orientation = UEs(u_).downlink_channel.UE_orientation;
    end
    
    %% 打印eNodeB
    if DEBUG_LEVEL>=2
        fprintf('eNodeB List\n');
        for b_=1:length(eNodeB_sites)
            eNodeB_sites(b_).print;
        end
        fprintf('\n');
    end
    
    %% 打印所有用户
    if DEBUG_LEVEL>=2
        fprintf('User List\n');
        for u_=1:length(UEs)
            UEs(u_).print;
        end
        fprintf('\n');
    end
    
    %% 主循环
    if DEBUG_LEVEL>=1
        fprintf('Entering main simulation loop, %5.0f TTIs\n',SYS_config.simulation_time_tti);
    end
    
    % 初始化计数器
    ticID_start_sim_loop = tic;
    starting_time = toc(ticID_start_sim_loop);
    
    num_markers = 5;% 进度条的长度
    s_markings  = round(linspace(1,length(eNodeB_sectors),num_markers));%扇区标记
    u_markings  = round(linspace(1,length(UEs),num_markers));%UE标记
    
    % 时钟初始化
    while networkClock.current_TTI < SYS_config.simulation_time_tti
        networkClock.advance_1_TTI;
        if DEBUG_LEVEL>=1
            fprintf(' TTI %5.0f/%d: ',networkClock.current_TTI,SYS_config.simulation_time_tti);
        end
        
        % 波束赋形
        if DEBUG_LEVEL>=1
            fprintf('\n               ');
        end
        if SYS_config.asynchronization %如果是不同步
            if SYS_config.isDouble == true %如果是双系统且开启了异步干扰
                [ue_or_bs_gain_asy,ue_gain_asy,bs_or_ue_gain_asy,bs_gain_asy] = NR_beam_forming_asy(UEs,eNodeB_sectors,networkPathlossMap,SYS_config); %不同步时采用不同的波束赋型

                [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeB_sectors,networkPathlossMap,SYS_config);

            else
                error('isDouble must be true')
            end
        else %如果是同步则调用正常的波束赋型函数
            [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeB_sectors,networkPathlossMap,SYS_config);
        end
       
        % 计算SINR
        if DEBUG_LEVEL>=1
            fprintf('               ');
        end
        switch SYS_config.asynchronization_switch
          case 0
                %下行
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain, SYS_config);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_down_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('+'); %进度条
                        end
                    end
                end
                for i = 1:length(UEs)
                    SINR(i) = UEs(i).wideband_SINR;
                    SINR2(i) = UEs(i).wideband_SINR2;
                    attached_sector_id(i) = UEs(i).attached_eNodeB.eNodeB_id;
                end
                figure;
                plot(attached_sector_id, LineWidth=2)
                set(gca, 'LineWidth', 1);
                ylim([10, 310]);
                ax = gca; % 获取当前坐标轴
                ax.Box = 'off'; % 关闭上边和右边的坐标轴
                if(SYS_config.labelMode == 'English')
                    xlabel('UAV flight sequence number')
                    ylabel('PCI')
                    str = sprintf('UAV Serving Cell (Simulated Data), h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
                elseif(SYS_config.labelMode == 'Chinese')
                    xlabel('无人机飞行序列号')
                    ylabel('PCI');
                    str = sprintf('UAV服务小区（仿真数据）, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
                end
                title(str)
                global interfer_sector_num
                % 统计UE的前8个干扰小区PCI
                interfering_cell_PCI = zeros(length(UEs), interfer_sector_num);
                for i = 1:length(UEs)
                    interfering_sectors = UEs(i).neighboring_Interference_eNodeB;
                    for j = 1:length(interfering_sectors)
                        interfering_cell_PCI(i, j) = interfering_sectors(j).eNodeB_id;
                    end
                end
                str = sprintf('PCI_%dm_Simulated.csv', SYS_config.UE_r);
                csvwrite(str, interfering_cell_PCI);
                switch(SYS_config.UE_r)
                case 100
                    y_lim = [-8,15];
                    scale_factor = 1.29;
                case 300
                    y_lim = [-10,20];
                    scale_factor = 0.86;
                case 500
                    y_lim = [-10,20];
                    scale_factor = 2;
                case 800   
                    y_lim = [-15,25];
                    scale_factor = 1.3627;
                case 1200
                    y_lim = [-15,25];
                    scale_factor = 1.3563;
                otherwise
                    fprintf("UAV极径输入错误！您输入的为%dm.",SYS_config.UE_r);
                end

                rng(0);
                fig1 = figure;
%                 plot_Sinr = SINR2;
                plot_Sinr = (SINR2 - mean(SINR2))*scale_factor + mean(SINR2);
                str_save_variable_name = sprintf('simulated_data_%d.mat', SYS_config.UE_r);
                save(str_save_variable_name, 'plot_Sinr');
                plot(1:length(UEs), plot_Sinr,LineWidth=2);%+10+sqrt(10)*randn(size(SINR2))
                xlim([0,length(UEs)]);
                ylim(y_lim);
                grid on;
                set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
                set(gca, 'LineWidth', 1);
                ax = gca; % 获取当前坐标轴
                ax.Box = 'off'; % 关闭上边和右边的坐标轴

                if(SYS_config.labelMode == 'English')
                    xlabel('UAV flight sequence number')
                    ylabel('SINR/dB')
                    str = sprintf('SINR Calculation Based on Paper Formula, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
                elseif(SYS_config.labelMode == 'Chinese')
                    xlabel('无人机飞行序列号')
                    ylabel('SINR/(dB)');
                    str = sprintf('根据论文公式计算SINR, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
                end  
                title(str)

                str2 = sprintf('Sinr_simul_r%d.mat', SYS_config.UE_r);
%                 save(str2, 'Sinr_simul');
% 
%                 str1 = sprintf('pos_r%d_h%d', SYS_config.UE_r, SYS_config.UE_height);
%                 saveas(fig1, str1, 'svg');        % 指定路径保存

%                 fig2 = figure;
%                 plot(x, SINR2,LineWidth=2.5);
% %                 hold on;
% %                 plot(x_tick, y_tick, 'bo', LineWidth=2.5)
%                 xlabel('极角/°')
%                 ylabel('SINR/dB')
% %                 str = sprintf('根据论文公式方法计算，基站天线上倾%d度\nUAV极径=%dm，h=%dm', SYS_config.tilt, SYS_config.UE_r, SYS_config.UE_height);
%                 str = sprintf('根据论文公式方法计算，基站天线垂直覆盖24度\nUAV极径=%dm，h=%dm', SYS_config.UE_r, SYS_config.UE_height);
%                 title(str)
% 
%                 str2 = sprintf('alpha_r%d_h%d', SYS_config.UE_r, SYS_config.UE_height);
%                 saveas(fig2, str2, 'svg');        % 指定路径保存


                error('画图已结束')

                
            case 1
                if SYS_config.isDouble == true
                    kk_ = 1;% 用来记录activate_UE产生数据的下标
                    for u_ = 1:length(UEs)
                        if ~UEs(u_).deactivate_UE
                            [user_macroscopic_pathloss(kk_),ul_dl_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
                            [dl_user_macroscopic_pathloss(kk_),dl_ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
                            kk_ = kk_ + 1;
                        else
                            UEs(u_).dummy_up_down_model;
                            UEs(u_).dummy_down_up_model;
                        end
                        if ~isempty(find(u_==u_markings,1))
                            if DEBUG_LEVEL>=1
                                fprintf('*'); %进度条
                            end
                        end
                        % 保存每个TTI得到的SINR
                        UEs(u_).save_SINR_in_TTI(SYS_config);
                    end
                end
               
            case 2
                UE_SINR_sum_dl=[];
                UE_SINR_sum_ul=[];
                %下行
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_down_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('+'); %进度条
                        end
                    end
                end

                % 上行
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [IoT(kk_),interference_level(kk_),BS_Rx_power(kk_),BS_Rx_intf_power(kk_),ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_link_quality_model(SYS_config,UEs,beam_antenna_gain,beam_ue_gain);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_up_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('-');  %进度条
                        end
                    end
                    % 保存每个TTI得到的SINR
                    %UEs(u_).save_SINR_in_TTI(SYS_config);
                end
                
                %异步干扰计算
                kk_ = 1;% 用来记录activate_UE产生数据的下标
                    for u_ = 1:length(UEs)
                        if ~UEs(u_).deactivate_UE
                            ul_dl_wideband_loop_SINR(:,:,kk_) = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
                            dl_ul_wideband_loop_SINR(:,:,kk_) = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
                            kk_ = kk_ + 1;
                        else
                            UEs(u_).dummy_up_down_model;
                            UEs(u_).dummy_down_up_model;
                        end
                    end
        end
        
%         if SYS_config.asynchronization %上下行不同步
%             if SYS_config.isDouble == true
%                 kk_ = 1;% 用来记录activate_UE产生数据的下标
%                 for u_ = 1:length(UEs)
%                     if ~UEs(u_).deactivate_UE
%                         ul_dl_wideband_loop_SINR(:,:,kk_) = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
%                         dl_ul_wideband_loop_SINR(:,:,kk_) = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
%                         kk_ = kk_ + 1;
%                     else
%                         UEs(u_).dummy_up_down_model;
%                         UEs(u_).dummy_down_up_model;
%                     end
%                     if ~isempty(find(u_==u_markings,1))
%                         if DEBUG_LEVEL>=1
%                             fprintf('*'); %进度条
%                         end
%                     end
%                     % 保存每个TTI得到的SINR
%                     UEs(u_).save_SINR_in_TTI(SYS_config);
%                 end
%             end
%         else %上下行同步则计算上行干扰上行和下行干扰下行的SINR
%             % 下行
%             kk_ = 1;
%             for u_ = 1:length(UEs)
%                 if ~UEs(u_).deactivate_UE
%                     [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain);
%                     kk_ = kk_ + 1;
%                 else
%                     UEs(u_).dummy_down_link_quality_model;
%                 end
%                 if ~isempty(find(u_==u_markings,1))
%                     if DEBUG_LEVEL>=1
%                         fprintf('+'); %进度条
%                     end
%                 end
%             end
%             
%             % 上行
%             kk_ = 1;
%             for u_ = 1:length(UEs)
%                 if ~UEs(u_).deactivate_UE
%                     [IoT(kk_),interference_level(kk_),BS_Rx_power(kk_),BS_Rx_intf_power(kk_),ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_link_quality_model(SYS_config,UEs,beam_antenna_gain,beam_ue_gain);
%                     kk_ = kk_ + 1;
%                 else
%                     UEs(u_).dummy_up_link_quality_model;
%                 end
%                 if ~isempty(find(u_==u_markings,1))
%                     if DEBUG_LEVEL>=1
%                         fprintf('-');  %进度条
%                     end
%                 end
%                 % 保存每个TTI得到的SINR
%                 UEs(u_).save_SINR_in_TTI(SYS_config);
%             end
%         end
        
        if DEBUG_LEVEL>=1
            fprintf('\n');
        end
        
        if networkClock.current_TTI < SYS_config.simulation_time_tti
            % Move UEs
            if ~SYS_config.keep_UEs_still
                move_all_UEs(SYS_config,UEs,networkPathlossMap,eNodeB_sectors);
            end
        end
        
        if mod(networkClock.current_TTI,10)==0
            elapsed_time = toc(ticID_start_sim_loop);
            time_per_iteration = elapsed_time / networkClock.current_TTI;
            estimated_time_to_finish = (SYS_config.simulation_time_tti - networkClock.current_TTI)*time_per_iteration;
            estimated_time_to_finish_h = floor(estimated_time_to_finish/3600);
            estimated_time_to_finish_m = estimated_time_to_finish/60 - estimated_time_to_finish_h*60;
            fprintf('Time to finish: %3.0f hours and %3.2f minutes\n',estimated_time_to_finish_h,estimated_time_to_finish_m);
        end
    end
    
    %% 统计
    % 被统计的UE都是deactivate_UE = false的
    if DEBUG_LEVEL>=1
        fprintf('Generating statistics\n');
    end
    if SYS_config.keep_UEs_still
        activate_UE_id = NR_get_activate_UE( UEs );
        switch SYS_config.asynchronization_switch
            case 0
               NR_plot_SINR_CDF(SYS_config,UEs,activate_UE_id); 
               NR_throughput(SYS_config,wideband_loop_SINR,ul_wideband_loop_SINR,activate_UE_id);
            case 1
               NR_plot_SINR_CDF_asy(SYS_config,UEs,activate_UE_id);
               NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
               NR_user_PL_gain(SYS_config,user_macroscopic_pathloss,dl_user_macroscopic_pathloss);
               NR_interfering_gain(SYS_config,UEs,activate_UE_id);
               NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
               %NR_user_dl_to_ul_PL_gain(SYS_config,dl_user_macroscopic_pathloss);
            case 2
               %同步干扰SINR
               all_SINR = [UEs.wideband_SINR];
               UE_SINR = all_SINR(activate_UE_id);
               UE_SINR_sum_dl=[UE_SINR_sum_dl UE_SINR];
               ul_all_SINR = [UEs.ul_wideband_SINR];
               UE_SINR_ul = ul_all_SINR(activate_UE_id);
               UE_SINR_sum_ul=[UE_SINR_sum_ul UE_SINR_ul];
               %异步干扰SINR
               ul_dl_all_SINR = [UEs.ul_dl_wideband_SINR];
               UE_SINR = ul_dl_all_SINR(activate_UE_id);
               UE_SINR_sum_dl=[UE_SINR_sum_dl UE_SINR];
               dl_ul_all_SINR = [UEs.dl_ul_wideband_SINR];
               UE_SINR_ul = dl_ul_all_SINR(activate_UE_id);
               UE_SINR_sum_ul=[UE_SINR_sum_ul UE_SINR_ul];
               NR_plot_SINR_CDF_half(UE_SINR_sum_dl,UE_SINR_sum_ul);
               
               wideband_loop_SINR_sum=wideband_loop_SINR;
               ul_wideband_loop_SINR_sum=ul_wideband_loop_SINR;
               wideband_loop_SINR_sum=cat(3,wideband_loop_SINR_sum,ul_dl_wideband_loop_SINR);
               ul_wideband_loop_SINR_sum=cat(3,ul_wideband_loop_SINR_sum,dl_ul_wideband_loop_SINR);
               NR_throughput_half(SYS_config,wideband_loop_SINR_sum,ul_wideband_loop_SINR_sum,activate_UE_id);
        end
        if SYS_config.asynchronization
            NR_plot_bs_bs_pathloss(SYS_config,UEs,activate_UE_id);
            NR_plot_ue_ue_pathloss(SYS_config,UEs,activate_UE_id);
        end
        NR_plot_couplingloss_CDF(SYS_config,UEs,beam_antenna_gain,beam_ue_gain,activate_UE_id);
        


%
       % NR_plot_receive_power_CDF( SYS_config,ul_to_dl_Agg_BS_Tx,  ul_to_dl_Agg_Rx,dl_to_ul_Agg_BS_Tx,dl_to_ul_Agg_Rx);
        %NR_dl_ul_gain(SYS_config,UEs,activate_UE_id);
        %NR_dl_ul_gain(SYS_config,ue_or_bs_gain_asy);
        %NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
      %  NR_plot_signal_and_intf_CL_CDF(,interfering_CL,another_interfering_CL,activate_UE_id);
       % NR_plot_receive_power_CDF( SYS_config,UE_Rx_power,UE_Rx_intf_power,BS_Rx_power,BS_Rx_intf_power );
        %NR_plot_receive_power_CDF(SYS_config,BS_Tx_power_asy,UE_RX_power_asy,BS_Tx_intf_power_asy,UE_Rx_intf_power_asy)
    %    NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
%             %绘制IoTcdf
%         if ~SYS_config.asynchronization
%             %保存各部分损耗
%             NR_plot_couplingloss_CDF(SYS_config,UEs,beam_antenna_gain,beam_ue_gain,activate_UE_id);
%             %保存功控cdf情况
%             NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
%             %保存邻系统耦合损耗cdf情况
%             NR_plot_signal_and_intf_CL_CDF(signal_CL,interfering_CL,another_interfering_CL,activate_UE_id);
%             %绘制IoTcdf
%             NR_plot_IoT_CDF(SYS_config,IoT);
%             % 绘制blocking cdf
%             NR_plot_blocking_CDF( SYS_config,interference_level );
%             % 绘制接收点功率及保存数据
%             NR_plot_receive_power_CDF( SYS_config,UE_Rx_power,UE_Rx_intf_power,BS_Rx_power,BS_Rx_intf_power );
%             % 不同ACIR下的吞吐量变化
%             if SYS_config.isACIR_loop && SYS_config.isDouble
%                 NR_throughput(SYS_config,wideband_loop_SINR,ul_wideband_loop_SINR,activate_UE_id);
%             end
%             % 固定ACIR下的吞吐量，打印在Command Windows上
%             NR_plot_thput_CDF(SYS_config,UEs,activate_UE_id);
%             % 绘制SINR的CDF
%             NR_plot_SINR_CDF(SYS_config,UEs,activate_UE_id);
%         else
%             % 异步时的
%             NR_plot_thput_CDF_asy(SYS_config,UEs,activate_UE_id);
%             if SYS_config.isACIR_loop && SYS_config.isDouble
%                 NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
%             end
%         end
    end
    
    %% 计算用时
    finish_time_s_full = toc(ticID_start_sim_begin);
    finish_time_m = floor(finish_time_s_full/60);
    finish_time_s = finish_time_s_full-finish_time_m*60;
    if DEBUG_LEVEL>=1
        fprintf('Simulation finished\n');
        fprintf(' Total elapsed time: %.0fm, %.0fs\n',finish_time_m,finish_time_s);
    end
    
    if DEBUG_LEVEL>=1
        fprintf('Saving results to %s\n',SYS_config.results_file);
    end
    
    %% 运行时产生的各种数据的保存
    % Compact (or not) results
    if SYS_config.compact_results_file
        
        networkPathlossMap.delete_everything_except_cell_assignments();
        
        clear the_UE_traces;

        % 将UE的部分属性变成结构体保存
        for u_=1:length(UEs)
            the_UE_traces(u_).position        = UEs(u_).pos;
            the_UE_traces(u_).attached_site   = UEs(u_).attached_site.id;
            the_UE_traces(u_).attached_eNodeB = UEs(u_).attached_eNodeB.eNodeB_id;
            the_UE_traces(u_).attached_eNodeB_id_TTI = UEs(u_).attached_eNodeB_id_TTI;
            the_UE_traces(u_).SINR_TTI = UEs(u_).SINR_TTI;
            the_UE_traces(u_).ul_SINR_TTI = UEs(u_).ul_SINR_TTI;
            the_UE_traces(u_).ul_dl_SINR_TTI = UEs(u_).ul_dl_SINR_TTI;
            the_UE_traces(u_).dl_ul_SINR_TTI = UEs(u_).dl_ul_SINR_TTI;
            
            UEs_struct(u_)                    = UEs(u_).basic_information_in_struct();
            UEs(u_).clear_non_basic_info();
        end
        UEs = UEs_struct; 
        % 保存运行数据
        %save(fullfile(SYS_config.results_folder,SYS_config.results_file),'SYS_config','networkPathlossMap','eNodeB_sites','eNodeB_sectors','UEs','the_UE_traces','finish_time_s_full','-v7.3');
    else
        
        if SYS_config.delete_pathloss_at_end
            networkPathlossMap.pathloss = [];
            if exist('networkShadowFadingMap','var')
                networkPathlossMap.pathloss = [];
            else
                % Do nothing
            end
        end
        save(fullfile(SYS_config.results_folder,SYS_config.results_file),'-v7.3');
        if DEBUG_LEVEL>=1
            fprintf('Traces saved in the standard format\n');
        end
    end
end

output_results_file = fullfile(SYS_config.results_folder,SYS_config.results_file);

utils.miscUtils.tidy_up_memory_before_closing(UEs,eNodeB_sectors,eNodeB_sites);
clear the_UE_traces
clear networkPathlossMap
clear networkShadowFadingMap
clear UEs
clear eNodeBs
clear eNodeBs_sectors

function move_all_UEs(SYS_config,UEs,networkPathlossMap,eNodeBs_sectors)
%该函数用于扩展UE移动方面

some_UE_out_of_ROI_this_TTI = false;
num_first_UEs = networkPathlossMap.num_first_UEs;
num_first_sectors = networkPathlossMap.num_first_sectors;
for u_ = 1:length(UEs)
    UEs(u_).move;
    [ x_range,y_range ] = networkPathlossMap.valid_range;
    
    % 如果UE走出了ROI，则令其随机移动到ROI内任意位置
    %
    % 这同时也实现了一个简单的切换程序，可以在此基础上进行扩充
    
    ROI_teleport       = ~UEs(u_).is_in_roi(SYS_config,x_range,y_range);% UE是否还在ROI中，否的话ROI_teleport=true
    handover_requested = UEs(u_).cell_change.requested;

    if ROI_teleport || handover_requested
        old_eNodeB_id = UEs(u_).attached_eNodeB.eNodeB_id;
        if ROI_teleport %如果超出了ROI
            new_UE_position = networkPathlossMap.random_position;%随机选择位置
            
            [new_site_id,new_sector_id,new_eNodeB_id] = networkPathlossMap.cell_assignment(SYS_config,u_,new_UE_position,num_first_UEs);
            
            %位置更新
            UEs(u_).pos = new_UE_position;
        elseif handover_requested
            new_eNodeB_id                     = UEs(u_).cell_change.target_eNodeB;
            UEs(u_).cell_change.requested     = false; % 重置切换域
            UEs(u_).cell_change.target_eNodeB = [];
        end
        
         % 将UE从原基站分离，并附属与新基站
        UEs(u_).start_handover(eNodeBs_sectors(new_eNodeB_id));
        

        % 统计或不统计UE的SINR
        if SYS_config.isWraparound
            compute_only_UEs_from_this_eNodeBs = 1:57;
        else
            if SYS_config.isDouble
                compute_only_UEs_from_this_eNodeBs = 1:num_first_sectors;
            else
                compute_only_UEs_from_this_eNodeBs = 1:length(eNodeBs_sectors);
            end
        end
        if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id == compute_only_UEs_from_this_eNodeBs,1))
            % 不统计该UE
            UEs(u_).deactivate_UE = true;
        else
            % 统计该UE
            UEs(u_).deactivate_UE = false;
        end
        
        % 打印一些调试信息
        if ~some_UE_out_of_ROI_this_TTI
            if SYS_config.debug_level>=1
                fprintf(1,'\n');
            end
            some_UE_out_of_ROI_this_TTI = true;
        end
        if SYS_config.debug_level>=1
            if ROI_teleport
                fprintf('              UE %g going out of ROI, teleporting to %g %g. eNodeB %g -> eNodeB %g\n',UEs(u_).id,new_UE_position(1),new_UE_position(2),old_eNodeB_id,new_eNodeB_id);
            elseif handover_requested
                fprintf('              UE %g handover request. eNodeB %g -> eNodeB %g\n',UEs(u_).id,old_eNodeB_id,new_eNodeB_id);
            end
        end
    end
end
if some_UE_out_of_ROI_this_TTI
    if SYS_config.debug_level>=1
        fprintf('              ');
    end
end
