classdef resourceBlockGrid < handle
    % RB类，主要作用是给每个RB分配功率
    % 现在只是用到将其作为一种定义使用带宽的方法，
    % 但是考虑到将来不是给每个UE全部的RB或者是上行ACIR模型像LTE里一样的话，
    % 这样建一个类是有必要的
    
    properties
        % BS发UE收的功率分配，单位w
        power_allocation
        % UE发BS收的功率分配，单位w
        ue_power_allocation
        % 导频的
        power_allocation_signaling
        ue_power_allocation_signaling
        % RB数
        n_RB
    end
    
    methods
        
        % 构造函数
        function obj = resourceBlockGrid(n_RB)
            obj.power_allocation = zeros(n_RB,1); 
            obj.power_allocation_signaling = zeros(n_RB,1);
            obj.ue_power_allocation = zeros(n_RB,1);
            obj.ue_power_allocation_signaling = zeros(n_RB,1);
            obj.n_RB = n_RB;
        end
        
        % 每个RB分配到的功率是一样的
        function set_homogeneous_power_allocation(obj,SYS_config,power_in_watts_data,power_in_watts_signaling)
            % power_in_watts_data = eNodeBs(b_).sectors(s_).max_power
            % power_in_watts_signaling = eNodeBs(b_).sectors(s_).signaling_power
            
            % 基站侧
            power_per_RB_data                 = power_in_watts_data / obj.n_RB;%功率平均分配给每个RB
            power_per_RB_signaling            = power_in_watts_signaling / obj.n_RB;
            obj.power_allocation(:)           = power_per_RB_data;
            obj.power_allocation_signaling(:) = power_per_RB_signaling;
            
            UE_tx_power = SYS_config.UE_tx_power;
            ue_max_power = UE_tx_power * (1-SYS_config.signaling_ratio);
            ue_signaling_power = UE_tx_power * SYS_config.signaling_ratio;
            
            % UE侧
            % 实际上基站与UE只是功率不同
            ue_power_per_RB_data = ue_max_power / obj.n_RB;%功率平均分配给每个RB
            ue_power_per_RB_signaling = ue_signaling_power / obj.n_RB;
            obj.ue_power_allocation(:)           = ue_power_per_RB_data;
            obj.ue_power_allocation_signaling(:) = ue_power_per_RB_signaling;
        end

        % Prints some info about this object
        function print(obj)
            fprintf('n_RB=%d\n',obj.n_RB);
        end
        
    end
end
