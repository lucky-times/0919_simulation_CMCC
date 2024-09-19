classdef macroscopicPathlossMap < handle
    % 路损图谱
    % 主要储存路损等各种参数，也包括了调用路损的函数
    % 同时也包含了记录小区分布的变量以及调用改变量的函数 
    
    properties
        % 路损包含天线增益
        pathloss
        % BS端天线增益
        BS_antenna_gain
         %UE天线增益
        UE_antenna_gain
        % 900直接计算得到的路损
        pathloss_900
        %NLOS判断
        LOS_compare
        %UE的高度
        H_height
        %UE的旋转角
        H_orientation
        % 穿透损耗
        penetration_loss
        % 距离
        distance_matrix
        micro_centre_distance_matrix
        % UE位置
        UE_pos_vector
        % 异构基站数
        num_hetnet_sites
        % 异构基站扇区数
        num_hetnet_sectors
         % 第一个系统基站数
        num_first_sites
        % 第一个系统扇区数
        num_first_sectors
        % 第一个系统UE数
        num_first_UEs
        
        % 第一个系统的sector_assignment
        % sector_assignment是一个二维数组，举个例子sector_assignment（3,5）=11，
        % 位于3行5列的那个像素点他的服务基站id为11
        sector_assignment
        
        % 第二个系统的sector_assignment
        sector_assignment_double
        
        % 第一个系统的sector_sizes，表示每个小区服务的像素点数
        sector_sizes
        
        % 第二个系统的sector_sizes
        sector_sizes_double
        
        % 地图分辨率
        data_res
        
        % 地图的XY轴
        roi_x
        roi_y
        
        % 描述这个路损模型的名字
        name
        
        sector_idx_mapping  % sector_idx_mapping(s_idx)  = [b_ s_]
        site_sector_mapping % site_sector_mapping(b_,s_) = s_idx

        % 扇区中心，用于GUI中标记
        sector_centers
        
        % 第二个系统的sector_centers
        sector_centers_double
       

    end

    methods
        function print(obj)
            fprintf('macroscopicPathlossMap\n');
            fprintf('Data resolution: %d meters/pixel\n',obj.data_res);
            fprintf('ROI: x: %d,%d y:%d,%d\n',obj.roi_x(1),obj.roi_x(2),obj.roi_y(1),obj.roi_y(2));
        end
        
        %% 下面一系列函数的作用是根据UE所在的位置返回保存在该类中的路损pathloss、UE高度H_height、UE旋转角度H_orientation等值
        % 方法大体相似：先得到UE的绝对位置，将其转换为像素点的位置（该像素点在矩阵中的索引），然后根据索引在矩阵调用
        
        %返回天线增益
        function BS_antenna_gain= get_bs_antenna_gain(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
%             if point_outside_lower_bound || point_outside_upper_bound
%                 BS_antenna_gain = NaN;
%             else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
               BS_antenna_gain = obj.BS_antenna_gain(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
%             end
        end
        
 %返回天线增益
        function UE_antenna_gain= get_ue_antenna_gain(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                UE_antenna_gain = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
               UE_antenna_gain = obj.UE_antenna_gain(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        
        % 返回路损
        function pathloss = get_pathloss(obj,pos,s_,b_)
            s_idx = obj.site_sector_mapping(b_,s_);
            s_idx = s_idx(:);
            s_idx = s_idx(s_idx~=0);
            pathloss = obj.get_pathloss_eNodeB(pos,s_idx);
        end
        % 返回某个点到某基站的路损
        function pathloss = get_pathloss_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                pathloss = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
                pathloss = obj.pathloss(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % 返回UE高度
        function height = get_UE_height(obj,pos)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                height = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
                height = obj.H_height(pixel_coord(:,2),pixel_coord(:,1));
            end
        end
        
        % 返回UE旋转角度
         function orientation = get_UE_orientation(obj,pos)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                orientation = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
                orientation = obj.H_orientation(pixel_coord(:,2),pixel_coord(:,1));
            end
         end
        
        % 返回没有天线增益的路损
        function pathloss_900 = get_pathloss_900_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                pathloss_900 = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
                pathloss_900 = obj.pathloss_900(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % 返回穿透损耗
        function penetration_loss = get_penetration_loss_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                penetration_loss = 0;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%就是将pos_to_pixel
                
                penetration_loss = obj.penetration_loss(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % Plots the pathloss of a given eNodeB sector
        function plot_pathloss(obj,b_,s_)
            figure;
            s_idx = obj.site_sector_mapping(b_,s_);
            imagesc(obj.pathloss(:,:,s_idx));
        end
        
        % Range of positions in which there are valid pathloss values
        function [x_range y_range] = valid_range(obj)
            x_range = [ obj.roi_x(1) obj.roi_x(2) ];
            y_range = [ obj.roi_y(1) obj.roi_y(2) ];
        end
        
        % Returns the coordinate origin for this pathloss map
        function pos = coordinate_origin(obj)
            pos = [ obj.roi_x(1) obj.roi_y(1) ];
        end
        
        % 给LTE系统加上MCL
        function apply_MCL(obj,sectors,minimum_coupling_loss)
            obj.pathloss(:,:,sectors) = max(obj.pathloss(:,:,sectors),minimum_coupling_loss);
        end
        
        % 知道UE的位置，返回该UE所属的扇区、基站
        function [ b_ s_ s_idx] = cell_assignment(obj,SYS_config,UE_id,pos,num_first_UEs,varargin)
            if length(varargin)<1
                x_ = pos(1);
                y_ = pos(2);
            else
                x_ = pos(1);
                y_ = varargin{1};
            end
            if max([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]) >= 1
                eNodeB_id = NaN;
                sector_num = NaN;
            elseif max([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]) >= 1
                eNodeB_id = NaN;
                sector_num = NaN;
            else
                pixel_coord = NR_common_pos_to_pixel([x_ y_],obj.coordinate_origin,obj.data_res);
                if SYS_config.isDouble
                    if UE_id<=num_first_UEs
                        s_idx = obj.sector_assignment(pixel_coord(2),pixel_coord(1));
                    else
                        s_idx = obj.sector_assignment_double(pixel_coord(2),pixel_coord(1));
                    end
                else
                    s_idx = obj.sector_assignment(pixel_coord(2),pixel_coord(1));%找到该UE所在位置对应的扇区，然后再找到对应的enb
                end
                
                b_s_        = obj.sector_idx_mapping(s_idx,:);
                b_          = b_s_(1);
                s_          = b_s_(2);
            end
        end
        
        % 返回基站数和扇区数
        function [ eNodeBs sectors_per_eNodeB ] = size(obj)
            eNodeBs            = size(obj.pathloss,4);
            sectors_per_eNodeB = size(obj.pathloss,3);
        end
        
        % 当UE运动移出了地图，重新给它分配一个位置
        function position = random_position(obj)
            [row,col] = find(obj.sector_assignment > 0);
            index = NR_randperm(length(col),1);
            position_pix = [col(index),row(index)];
            position = NR_common_pixel_to_pos(position_pix,obj.coordinate_origin,obj.data_res);
        end
        
        % Deletes all information except for the cell assignment maps
        function delete_everything_except_cell_assignments(obj)
            obj.pathloss            = [];
            obj.name                = [];
            obj.sector_idx_mapping  = [];
            obj.site_sector_mapping = [];
            obj.BS_antenna_gain     = [];
            obj.UE_antenna_gain     = [];
            obj.pathloss_900        = [];
            obj.LOS_compare         = [];
            obj.H_height            = [];
            obj.H_orientation       = [];
            obj.penetration_loss    = [];
            obj.distance_matrix     = [];
            obj.micro_centre_distance_matrix  = [];
            obj.UE_pos_vector       = [];
        end
    end
end