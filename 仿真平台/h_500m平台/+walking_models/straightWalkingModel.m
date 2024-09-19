classdef straightWalkingModel < walking_models.walkingModel
    % UE直线行走模型
    
    properties
        % 方向
        direction
        % 由UE_speed*TTI_length算出路程
        s
    end
    methods
        % 构造函数
        function obj = straightWalkingModel(s,varargin)
            % 如果没有指定方向则随机选择一个方向
            if length(varargin)<1
                direction = floor(random('unif',0,359));
            else
                direction = varargin{1};
            end
            obj.direction = direction;
            obj.s = s;
        end
        
        % 根据现在的位置算出下一个TTI的位置
        function new_pos = move(obj,current_pos)
            mov_vector = obj.s * [ cosd(obj.direction) sind(obj.direction) ];
            new_pos = current_pos + mov_vector;
        end
        
        % 根据现在的位置算出上一个TTI的位置
        function old_pos = move_back(obj,current_pos)
            mov_vector = obj.s * [ cosd(obj.direction) sind(obj.direction) ];
            old_pos = current_pos - mov_vector;
        end
        
        % Print some info
        function print(obj)
            fprintf('  direction: %d? %d m/s*TTI\n',obj.direction,obj.s);
        end
    end
end
