classdef eNodeB_sector < handle
    % 小区实体

    properties
        id                           % 小区在eNodeB中相对索引
        eNodeB_id                    % 本sector在所有sector中的PCI
        global_id                    % 本sector在所有sector中的id，不是本sector所属enb的id
        parent_eNodeB                % 小区所属的eNodeB
        azimuth                      % 小区天线扇瓣的旋转角
        antenna                      % 小区天线
        attached_UEs = 0;            % 小区连接的用户数
        max_power                    % 小区最大的发射功率
        signaling_power              % Power dedicated to signaling. Counted as always in use，为0
        scheduler                    % 小区的调度器
        RB_grid                      % 小区的资源块分配表
        nTX                          % 小区的（虚拟）天线数
        in_interf_eNodeB_sectors             % 第一个系统的干扰扇区
        ad_interf_eNodeB_sectors                % 第二个系统的干扰扇区
        
        always_on = true;            % 是否小区常开（无论是否有用户接入）
        
        last_received_feedback       % 最近接收信号反馈
        attached_UEs_vector          % 接入UE的列表
        
        feedback_trace               % 接入UE的列表
        sector_trace                 % 接入UE的列表
        
        zero_delay_feedback = false; % 接入UE的列表

        macroscopic_pathloss_model   % 小区使用的路损模型
        

        unquantized_CQI_feedback = false;
        
        % extra params that may not be always used (mainly extra information from the network planning tool)
        transmitter     % This points to the pathloss file that will actually be used
        frequency_band  % 小区使用的频带
        antenna_name    % 小区的天线名
        antenna_type    % 小区使用的天线类型
        tx_height       % 小区的天线高度
        electrical_downtilt % 小区天线的电气下倾角
        mechanical_downtilt % 小区天线的机械下倾角
    end

    methods
        
        function print(obj)
            fprintf(' Sector %d: ',obj.id);
            fprintf('%s %ddB %d\n',obj.antenna.antenna_type,obj.antenna.mean_antenna_gain,obj.azimuth);
            fprintf('  ');
            obj.scheduler.print;
            fprintf('  UEs: ');
            if ~isempty(obj.users)
                current_last_node = obj.users(1);
                while ~isempty(current_last_node.Next)
                    % Do something
                    fprintf('%d ',current_last_node.Data.id);
                    current_last_node = current_last_node.Next;
                end
                % Do something for the last node
                fprintf('%d ',current_last_node.Data.id);
            end
            fprintf('\n');
            fprintf('  '); obj.RB_grid.print;
        end
        
        function clear(obj)
            obj.parent_eNodeB              = [];
            obj.antenna                    = [];
            obj.attached_UEs_vector        = [];
            obj.scheduler                  = [];
            obj.RB_grid                    = [];
            obj.in_interf_eNodeB_sectors           = [];
            obj.ad_interf_eNodeB_sectors = [];
            obj.last_received_feedback     = [];
            obj.feedback_trace             = [];
            obj.sector_trace               = [];
            obj.macroscopic_pathloss_model = [];
            obj.transmitter                = [];
            obj.frequency_band             = [];
        end
        
        %% 将指定用户接入到小区中
        function attachUser(obj,user)
            if isempty(obj.attached_UEs_vector)
 
                obj.attached_UEs_vector  = user;
                obj.attached_UEs         = obj.attached_UEs + 1;
            else

                current_UEs = [obj.attached_UEs_vector.id];
                if ~sum(current_UEs==user.id)
                    obj.attached_UEs_vector = [obj.attached_UEs_vector user];
                    obj.attached_UEs        = obj.attached_UEs + 1;
                end
            end            
            % UE在这里attach上了
            user.attached_sector_idx = obj.id;
            user.attached_site       = obj.parent_eNodeB;
            user.attached_eNodeB     = obj;
        end
        %% 观察UE的干扰小区
        function neighboring_Interfer(obj, user)
            if isempty(user.neighboring_Interference_eNodeB)
                user.neighboring_Interference_eNodeB = obj;
            else
                user.neighboring_Interference_eNodeB = [user.neighboring_Interference_eNodeB obj];
            end
        end

        %% 将指定用户从小区中分离
        function deattachUser(obj,user)
            if ~isempty(obj.attached_UEs_vector)
                current_UEs  = [obj.attached_UEs_vector.id];
                UE_idx       = (current_UEs==user.id);
                UE_in_eNodeB = sum(UE_idx);
                
                if UE_in_eNodeB>0
                    obj.attached_UEs_vector = obj.attached_UEs_vector(~UE_idx);
                end
                obj.attached_UEs = obj.attached_UEs - UE_in_eNodeB;
            end
        end
        
         %% 判断某个用户是否接入该小区
        function is_attached = userIsAttached(obj,user)
            if ~isempty(obj.users)
                current_UEs  = [obj.attached_UEs_vector.id];
                UE_idx       = (current_UEs==user.id);
                UE_in_eNodeB = sum(UE_idx);
                is_attached  = logical(UE_in_eNodeB);
            else
                is_attached = false;
            end
        end       
        
        
        function clear_non_basic_info(obj)
            obj.antenna                    = [];
            obj.attached_UEs               = [];
            obj.attached_UEs_vector        = [];
            obj.scheduler                  = [];
            obj.RB_grid                    = [];
            obj.last_received_feedback     = [];
            obj.feedback_trace             = [];
            obj.sector_trace               = [];
            obj.macroscopic_pathloss_model = [];
        end
    end
end
