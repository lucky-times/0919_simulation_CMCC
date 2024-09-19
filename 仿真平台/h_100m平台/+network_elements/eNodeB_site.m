classdef eNodeB_site < handle
    % eNodeB实体（基站）
    
    properties
        id                  % eNodeB 的id
        pos                 % eNodeB的位置 (x,y)
        parent_centre_pos   % 主要存储umi扇区的中心
        sector_centre       % 主要存储uma扇区的中心
        femto_room          % 存储femto中房间大小
        O_InH;
        sectors             % eNodeB上的小区实体
        servering_height    % 基站服务高度
        in_interf_eNodeB_sites    % 第一个系统干扰基站
        ad_interf_eNodeB_sites       % 第二个系统干扰基站
        name
        altitude = 0;
        site_name
        
        site_type           
        clock               % 网络时钟
    end

    methods
        function print(obj)
            fprintf('eNodeB %d, position (x,y): (%3.2f,%3.2f), altitude: %3.0fm, %d attached UEs\n',obj.id,obj.pos(1),obj.pos(2),obj.altitude,obj.attached_UEs);
            fprintf('Sector: ');
            
            fprintf('\n');
            for s_=1:length(obj.sectors)
                obj.sectors(s_).print;
            end
        end
        
      
        function clear(obj)
            obj.sectors          = [];
            obj.in_interf_eNodeB_sites = [];
            obj.ad_interf_eNodeB_sites = [];
            obj.clock            = [];
        end
        
         %% 判断是否有用户连接到该eNodeB
        function is_attached = userIsAttached(obj,user)
            for s_ = 1:length(obj.sectors)
                if obj.sectors(1).userIsAttached(user);
                    is_attached = true;
                    return
                end
            end
            is_attached = false;
        end
        %% 得到该eNodeB连接的用户数 
        function number_of_atached_UEs = attached_UEs(obj)
            temp = zeros(1,length(obj.sectors));
            for s_ = 1:length(obj.sectors)
                temp(s_) = obj.sectors(s_).attached_UEs;
            end
            number_of_atached_UEs = sum(temp);
        end
    end
end