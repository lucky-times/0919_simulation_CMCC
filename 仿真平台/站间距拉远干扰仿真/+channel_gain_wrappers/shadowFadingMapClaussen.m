classdef shadowFadingMapClaussen < handle
    % 保存阴影衰落的类

    properties
        % 阴影衰落
        pathloss
        % 地图分辨率
        data_res
        % XY轴
        roi_x
        roi_y
    end
    
    methods
        % 构造函数
        function obj = shadowFadingMapClaussen(SYS_config,networkPathlossMap,eNodeB_sites)
            %% 和networkPathlossMap里的一样
            obj.roi_x       = networkPathlossMap.roi_x;
            obj.roi_y       = networkPathlossMap.roi_y;
            obj.data_res    = networkPathlossMap.data_res;
            obj.pathloss    = zeros(size(networkPathlossMap.pathloss));
            fprintf('No shadowFading!\n')
%             num_eNodeBs = length(eNodeB_sites);
%             [row,col,~] = size(networkPathlossMap.pathloss);
%             shadow_fading = zeros(row,col,num_eNodeBs);
%             LOS_compare = networkPathlossMap.LOS_compare;
%             r_eNodeBs = SYS_config.r_eNodeBs;
%             
%             %% 产生公共数据集
%             norm_public_LOS = normrnd(0, SYS_config.shadow_fading_sd_LOS,row,col); %产生LOS情况下的公共数据集
%             norm_public_NLOS = normrnd(0, SYS_config.shadow_fading_sd_NLOS,row,col); %产生NLOS情况下的公共数据
%             norm_public_NLOS_LTE = normrnd(0, SYS_config.shadow_fading_sd_NLOS_LTE,row,col);
%             switch SYS_config.scene_type
%                 case 'RMa'
%                     norm_public_LOS2 = normrnd(0, SYS_config.shadow_fading_sd_LOS2,row,col);
%                 case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
%                     norm_public_LOS_het = normrnd(0, SYS_config.shadow_fading_sd_LOS_het,row,col);
%                     norm_public_NLOS_het = normrnd(0, SYS_config.shadow_fading_sd_NLOS_het,row,col);
%             end
%             
%             %% 利用n_sites来判断不同基站产生不同的阴影衰落
%             if SYS_config.isDouble && SYS_config.shift_mode == 0
%                 n_sites = networkPathlossMap.num_first_sites;
%             else
%                 switch SYS_config.scene_type
%                     case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
%                         n_sites = networkPathlossMap.num_first_sites;
%                     otherwise
%                         n_sites = num_eNodeBs;
%                 end
%             end
%             % LTE与NR系统
%             if SYS_config.isDouble
%                 num_first_sites = networkPathlossMap.num_first_sites;
%                 first = class(eNodeB_sites(1).sectors(1).macroscopic_pathloss_model);
%                 second = class(eNodeB_sites(num_first_sites+1).sectors(1).macroscopic_pathloss_model);
%                 if ~strcmp(first,second)
%                     n_sites = num_eNodeBs;% 异路损模型一步算出
%                     isDiff = true;% 判断异步与否
%                 else
%                     isDiff = false;% 同路损模型不用n_sites赋值，因为是作为前面那些条件的补充
%                 end
%             end
%             
%             %% 单系统、双系统不共址、非异构、异路损模型一步直接算出所有基站的SF，双系统共址、异构算出第一个系统基站的SF
%             for b_ = 1:n_sites
%                 pathlossmodel = class(eNodeB_sites(b_).sectors(1).macroscopic_pathloss_model);
%                 if strcmp(pathlossmodel,'macroscopic_pathloss_models.TS38900PathlossModel')
%                     NR = true;% 用来判断是LTE还是5G
%                 else
%                     NR = false;
%                 end
%                 for r = 1:row
%                     for c = 1:col
%                         if LOS_compare(r,c,eNodeB_sites(b_).sectors(1).eNodeB_id) == 0
%                            shadow_fading(r,c,b_) = r_eNodeBs * norm_public_LOS(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_LOS);
%                         elseif LOS_compare(r,c,eNodeB_sites(b_).sectors(1).eNodeB_id) == 1
%                             if NR
%                                 shadow_fading(r,c,b_) = r_eNodeBs * norm_public_NLOS(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_NLOS);
%                             else
%                                 shadow_fading(r,c,b_) = r_eNodeBs * norm_public_NLOS_LTE(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_NLOS_LTE);
%                             end
%                         elseif LOS_compare(r,c,eNodeB_sites(b_).sectors(1).eNodeB_id) == -1 % RMa的LOS中有两种阴影衰落标准差
%                             shadow_fading(r,c,b_) = r_eNodeBs * norm_public_LOS2(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_LOS2);
%                         end
%                     end
%                 end
%             end
%             
%             %% 双系统共址或异构的第二个系统基站的SF
%             if SYS_config.isDouble && SYS_config.shift_mode == 0 % 共址基站,对UE采用相同的阴影衰落值
%                 if ~isDiff % 对应于双系统共址异路损模型的情形就不用这样赋值
%                     shadow_fading(:,:,n_sites+1:num_eNodeBs) = shadow_fading(:,:,1:n_sites);
%                 end
%             else % 异构的第二个系统基站的SF
%                 switch SYS_config.scene_type
%                     case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
%                         for b_ = n_sites+1:num_eNodeBs
%                             pathlossmodel = class(eNodeB_sites(b_).sectors(1).macroscopic_pathloss_model);
%                             if strcmp(pathlossmodel,'macroscopic_pathloss_models.TS38900PathlossModel')
%                                 NR = true;
%                             else
%                                 NR = false;
%                             end
%                             for r = 1:row
%                                 for c = 1:col
%                                     if LOS_compare(r,c,eNodeB_sites(b_).sectors(1).eNodeB_id) == 0
%                                         shadow_fading(r,c,b_) = r_eNodeBs * norm_public_LOS_het(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_LOS_het);
%                                     elseif LOS_compare(r,c,eNodeB_sites(b_).sectors(1).eNodeB_id) == 1
%                                         if NR
%                                             shadow_fading(r,c,b_) = r_eNodeBs * norm_public_NLOS_het(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_NLOS_het);
%                                         else
%                                             shadow_fading(r,c,b_) = r_eNodeBs * norm_public_NLOS_LTE(r,c) + sqrt(1-r_eNodeBs^2) * normrnd(0,SYS_config.shadow_fading_sd_NLOS_LTE);
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                 end
%             end
%             
%             obj.pathloss = shadow_fading;
            
        end

        function print(obj)
            fprintf('Data resolution: %d meters/pixel\n',obj.data_res);
            fprintf('ROI: x: %d,%d y:%d,%d\n',obj.roi_x(1),obj.roi_x(2),obj.roi_y(1),obj.roi_y(2));
        end
        
        % 返回某点的阴影衰落
        function pathloss = get_pathloss(obj,pos,b_)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                pathloss = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;
                pathloss = obj.pathloss(pixel_coord(:,2),pixel_coord(:,1),b_);
            end
        end
        
        % Plots the pathloss of a given eNodeB
        function plot_pathloss(obj,b_)
            figure;
            imagesc(obj.pathloss(:,:,b_));
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
        
        % Returns the number of eNodeBs that this pathloss map contains
        function [ eNodeB_sites ] = size(obj)
            eNodeB_sites            = size(obj.pathloss,3);
        end
    end
end