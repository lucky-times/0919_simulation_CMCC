classdef generalUeSpatialDistribution < handle
    % 该类是控制UE撒点的父类
    % 构造函数作用是把networkPathlossMap，eNodeB_sites，eNodeB_sectors传进去
    % 该类从以前的平台继承下来，以前的平台还有别的撒点方式，当我们的平台就只是用每个小区撒相同的UE数这种
    
    properties
        networkPathlossMap % 主要是调用其sector_assignment属性
        eNodeB_sites              % 基站
        eNodeB_sectors            % 扇区
    end
    
    methods
        function obj = generalUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors)
            obj.networkPathlossMap = networkPathlossMap;
            obj.eNodeB_sites              = eNodeB_sites;
            obj.eNodeB_sectors            = eNodeB_sectors;
        end
    end

    methods(Abstract)
        user_positions_pixels = generate_UE_positions(obj,varargin)
    end
    
end

