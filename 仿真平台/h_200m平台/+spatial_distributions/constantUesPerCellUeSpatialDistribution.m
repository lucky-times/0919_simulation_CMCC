classdef constantUesPerCellUeSpatialDistribution < spatial_distributions.generalUeSpatialDistribution
    % 用于撒点的一个类

    properties
        % 每个小区UE数
        UE_per_eNodeB
        % 异构下第二个系统每个小区UE数
        UE_per_eNodeB2
        % 双系统时第一个系统扇区数
        num_first_sectors
    end
    
    methods
        % Class constructor.
        function obj = constantUesPerCellUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors,SYS_config)
            % Fill in basic parameters (handled by the superclass constructor)
            obj = obj@spatial_distributions.generalUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors);%@用于派生类的构造函数，@用于调用基类成员函数
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                    obj.UE_per_eNodeB = SYS_config.UE_per_eNodeB;
                    obj.UE_per_eNodeB2 = SYS_config.UE_per_eNodeB2;
                    obj.num_first_sectors = length(eNodeB_sectors)-networkPathlossMap.num_hetnet_sectors;% 只是在双系统下
                otherwise
                    obj.UE_per_eNodeB = SYS_config.UE_per_eNodeB;
                    obj.UE_per_eNodeB2 = 0;
                    obj.num_first_sectors = length(eNodeB_sectors)/2;% 只是在双系统下
            end
        end
        
        function user_positions_pixels = generate_UE_positions(obj,SYS_config)
            networkPathlossMap = obj.networkPathlossMap;
            sector_surfaces    = networkPathlossMap.sector_sizes;
            sector_surfaces2 = networkPathlossMap.sector_sizes_double;
            eNodeB_sectors            = obj.eNodeB_sectors;%sector
            
            % 出现基站的服务范围为0的情况，出现这种情况一种可能是分辨率太粗，还有可能是前面哪里发生了错误
            % 当这种情况发生时会导致求SINR那里报错
            if sum(sector_surfaces(:)==0)>0
                warning('Some sector_surfaces are zero. No UEs will be generated there. Maybe a too big map resolution value?');
            end
            
            if SYS_config.isDouble
                if sum(sector_surfaces2(:)==0)>0
                    warning('Some sector_surfaces2 are zero. No UEs will be generated there. Maybe a too big map resolution value?');
                end
            end
            
            % users_sector是某个扇区要撒的UE数
            norm_sector_surface = ones(size(sector_surfaces));
            norm_sector_surface2 = ones(size(sector_surfaces2));
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}%没弄懂
                    users_sector = round(norm_sector_surface*obj.UE_per_eNodeB).*(sector_surfaces>0);
                    users_sector2 = round(norm_sector_surface2*obj.UE_per_eNodeB2).*(sector_surfaces2>0);
                otherwise
                    users_sector = round(norm_sector_surface*obj.UE_per_eNodeB).*(sector_surfaces>0);
                    users_sector2 = round(norm_sector_surface2*obj.UE_per_eNodeB).*(sector_surfaces2>0);
            end
            
            % sector_positions是扇区服务的所有像素点的位置，user_positions_pixels_per_cell是随机撒完点后像素点的位置
            sector_positions = cell(size(networkPathlossMap.sector_sizes));%cell是一种可以存储不同类型数据的数组
            user_positions_pixels_per_cell = cell(size(networkPathlossMap.sector_sizes));
            sector_positions2 = cell(size(networkPathlossMap.sector_sizes_double));
            user_positions_pixels_per_cell2 = cell(size(networkPathlossMap.sector_sizes_double));
            
            % Assign random positions to each UE
            if SYS_config.isDouble
                n_sector = obj.num_first_sectors;
            else
                n_sector = length(eNodeB_sectors);
            end
            for s_idx = 1:n_sector
                if users_sector(s_idx)~=0                   
                    [row,col] = find(networkPathlossMap.sector_assignment(:,:)==s_idx);
                    sector_positions{s_idx} = [col,row];%cell的大括号表示内容
                    % 在基站的服务范围内随机选出UE_per_eNodeB个点作为UE的位置
                    user_positions_pixels_per_cell{s_idx} = sector_positions{s_idx}(randi(size( sector_positions{s_idx},1)  ,  [1 users_sector(s_idx)]  ),:);
                else
                    sector_positions{s_idx}      = [];
                    user_positions_pixels_per_cell{s_idx} = [];
                end
            end
            if SYS_config.isDouble % 第二个系统调用不同的sector_assignment
                for s_idx = n_sector+1:length(eNodeB_sectors)
                    if users_sector2(s_idx-n_sector)~=0
                        [row,col] = find(networkPathlossMap.sector_assignment_double(:,:)==s_idx);
                        sector_positions2{s_idx-n_sector} = [col,row];%cell的大括号表示内容
                        % 在基站的服务范围内随机选出UE_per_eNodeB个点作为UE的位置
                        user_positions_pixels_per_cell2{s_idx-n_sector} = sector_positions2{s_idx-n_sector}(randi(size( sector_positions2{s_idx-n_sector},1)  ,  [1 users_sector2(s_idx-n_sector)]  ),:);
                    else
                        sector_positions2{s_idx-n_sector}      = [];
                        user_positions_pixels_per_cell2{s_idx-n_sector} = [];
                    end
                end
                user_positions_pixels_per_cell = [user_positions_pixels_per_cell user_positions_pixels_per_cell2];
            end
            
            % 下面这步是将上面得到的UE的位置整合成n_UE行，2列的形式，第一列是UE的x轴，第二列是UE的y轴
            user_positions_pixels_total_UEs = 0;
            for i_=1:length(user_positions_pixels_per_cell)
                if ~isempty(user_positions_pixels_per_cell{i_})
                    user_positions_pixels_total_UEs = user_positions_pixels_total_UEs + size(user_positions_pixels_per_cell{i_},1);
                end
            end
            user_positions_pixels = zeros(user_positions_pixels_total_UEs,2);
            l_=1;
            for i_=1:length(user_positions_pixels_per_cell)
                for j_=1:size(user_positions_pixels_per_cell{i_},1)
                    user_positions_pixels(l_,:) = user_positions_pixels_per_cell{i_}(j_,:);
                    l_=l_+1;
                end
            end
        end
    end
    
end

