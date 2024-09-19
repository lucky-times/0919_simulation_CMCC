function NR_calculate_cell_coverage(SYS_config,networkPathlossMap,eNodeB_sites,eNodeB_sectors,networkShadowFadingMap)
%基于系统设置初步计算小区容量
%输入参数：
%SYS_config：LTE参数配置系统
%networkPathlossMap：路损图谱
%eNodeB_sites：基站
%eNodeBs_sectors：扇区
%networkShadowFadingMap：阴影衰落映射矩阵图

if SYS_config.debug_level>=1
    fprintf('Calculating cell capacity:');
end

%% 求SNR
num_eNodeBs = length(eNodeB_sites);

shift_pathloss = networkPathlossMap.pathloss;
penetration_loss = networkPathlossMap.penetration_loss;
num_hetnet_sites = networkPathlossMap.num_hetnet_sites;
num_hetnet_sectors = networkPathlossMap.num_hetnet_sectors;
% UE端有天线增益，BS端没有天线增益，但是接入时仍然考虑BS天线增益
RX_powers_W = zeros(size(networkPathlossMap.pathloss));% 初始化接收功率矩阵
for b_ = 1:num_eNodeBs
    
    % 阴影衰落
    shadow_fading_current_eNodeB = 10.^(networkShadowFadingMap.pathloss(:,:,b_)/10);
    
    for s_ = 1:length(eNodeB_sites(b_).sectors)
        s_idx = eNodeB_sites(b_).sectors(s_).eNodeB_id;
        %计算接收功率
        RX_powers_W(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).max_power./10.^(shift_pathloss(:,:,s_idx)/10) ./ shadow_fading_current_eNodeB ./ 10.^(penetration_loss(:,:,s_idx)/10);
    end
end
%% 计算SNR
% 这里不管一个基站同时服务多少个UE，都认为基站对一个UE使用全部资源，因为这样不会影响求基站的服务区域
if SYS_config.debug_level>=1
    fprintf('calculating SNR-->');
end
SNR_linear_all = zeros(size(RX_powers_W));
thermal_noise_W = 10^(SYS_config.UE.thermal_noise_density/10) / 1000 * SYS_config.bandwidth * 10^(SYS_config.UE_receiver_noise_figure/10);

tot_eNodeBs = size(RX_powers_W,3);%第三维为小区个数
if SYS_config.isDouble
    num_first_sectors = networkPathlossMap.num_first_sectors;
    num_first_sites = networkPathlossMap.num_first_sites;
end
for s_=1:tot_eNodeBs
    SNR_linear_all(:,:,s_) = RX_powers_W(:,:,s_)/thermal_noise_W;
end
SNR_dB_all = 10*log10(SNR_linear_all);

% 对每个像素点接收到的SNR进行排序，选择接收到最大SNR时对应的基站作为服务基站进行接入
if SYS_config.isDouble
    [~,SNR_IX] = sort(SNR_dB_all(:,:,1:num_first_sectors),3);
    [~,SNR_IX2] = sort(SNR_dB_all(:,:,num_first_sectors+1:tot_eNodeBs),3);% +num_first_sectors是因为排序出来的索引是从1开始的
else
    [~,SNR_IX] = sort(SNR_dB_all,3);
end

SNR_assignment  = SNR_IX(:,:,end);%SNR_assignment返回每个像素点最大的SNR对应的扇区编号
if SYS_config.isDouble
    SNR_assignment2 = SNR_IX2(:,:,end)+num_first_sectors;% +num_first_sectors是因为排序出来的索引是从1开始的
else
    SNR_assignment2 = [];
end

% 3dB hand over
if SYS_config.hand_over
    if SYS_config.debug_level>=1
        fprintf('implement 3dB handover-->');
    end
    if SYS_config.isDouble
        % 第一个系统
        SNR_dB_all_tmp = SNR_dB_all(:,:,1:num_first_sectors);
        for s_ = 1:num_first_sectors
            [m,n] = find(SNR_assignment == s_);% 当前基站服务的像素点
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;% 这些像素点SNR-3
            [~,SNR_IX_tmp] = sort(SNR_dB_all_tmp,3);% -3后重新排序
            
            for i = 1:length(m)
                if SNR_IX_tmp(m(i),n(i),end) ~= s_ % 若这些像素点接收到最大SNR的基站发生了改变
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));% 找出比-3后SNR值大或等的基站
                    select = index(NR_randperm(length(index),1));% 从这些基站中随机选择一个作为服务基站
                    SNR_assignment(m(i),n(i)) = select;
                end
            end
        end
        % 第二个系统
        SNR_dB_all_tmp = SNR_dB_all(:,:,num_first_sectors+1:tot_eNodeBs);
        for s_ = 1:tot_eNodeBs-num_first_sectors
            [m,n] = find(SNR_assignment2 == (s_+num_first_sectors));
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;
            [~,SNR_IX2_tmp] = sort(SNR_dB_all_tmp,3);
            
            for i = 1:length(m)
                if SNR_IX2_tmp(m(i),n(i),end) ~= s_
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));
                    select = index(NR_randperm(length(index),1));
                    SNR_assignment2(m(i),n(i)) = select+num_first_sectors;
                end
            end
        end
        
    else
        SNR_dB_all_tmp = SNR_dB_all;
        for s_ = 1:tot_eNodeBs
            [m,n] = find(SNR_assignment == s_);% 当前基站服务的像素点
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;% 这些像素点SNR-3
            [~,SNR_IX_tmp] = sort(SNR_dB_all_tmp,3);% -3后重新排序
            
            for i = 1:length(m)
                if SNR_IX_tmp(m(i),n(i),end) ~= s_ % 若这些像素点接收到最大SNR的基站发生了改变
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));% 找出比-3后SNR值大或等的基站
                    select = index(NR_randperm(length(index),1));% 从这些基站中随机选择一个作为服务基站
                    SNR_assignment(m(i),n(i)) = select;
                end
            end
        end
        SNR_assignment2 = [];
    end
   
end

% 上面得到的SNR_assignment意思：假如SNR_assignment=(1,1,2,3,3,3),代表（1,1）这个地方由基站1服务，（1,5）这个地方由基站3服务

% 为了将UE撒在提案上的拓扑范围内，修改SINR_assignment
% 大体思路是：先将原来的SINR_assignment保存好后赋-1，-1代表不能在该处撒点
% 然后找出可以撒点的地方，对应处赋0
% 最后在SINR_assignment等于0的地方用原来的值替换
if SYS_config.debug_level>=1
    fprintf('restrict cell coverage\n');
end
switch SYS_config.scene_type
    case {'UMA','RMa'}
        ISD = SYS_config.ISD;
        SNR_assignment_tmp = SNR_assignment;% 先保存原来的SNR_assignment
        SNR_assignment(:,:) = -1;% 将SNR_assignment赋成-1，-1代表不能在此处撒点的意思
        distance_matrix = networkPathlossMap.distance_matrix;
        if SYS_config.isDouble
            n_enb = num_eNodeBs/2;
        else
            n_enb = num_eNodeBs;
        end
        for b_ = 1:n_enb
            [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);%UMA撒点的位置不能在基站的35m范围内
            for i = 1:length(m)
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                    if isinhex
                       SNR_assignment(m(i),n(i)) = 0; % 能够撒点的地方赋0
                    end
                end
            end
        end
        [m,n] = find(SNR_assignment(:,:)==0);
        for i = 1:length(m)
            SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));% 将可以撒点的地方用原来的值替换
        end
        if SYS_config.isDouble
            SNR_assignment_tmp2 = SNR_assignment2;
            SNR_assignment2(:,:) = -1;
            for b_ = n_enb+1:num_eNodeBs
                [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
                for i = 1:length(m)
                    for s_ = 1:length(eNodeB_sites(b_).sectors)
                        isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                        if isinhex
                            SNR_assignment2(m(i),n(i)) = 0;
                        end
                    end
                end
            end
            [m,n] = find(SNR_assignment2(:,:)==0);
            for i = 1:length(m)
                SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
            end
        end
    case 'UMI'
        if ~SYS_config.isManhattan
            UMi_r = SYS_config.UMi_r;
            SNR_assignment_tmp = SNR_assignment;
            SNR_assignment(:,:) = -1;
            distance_matrix = networkPathlossMap.distance_matrix;
            micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
            if SYS_config.isDouble
                n_enb = num_eNodeBs/2;
            else
                n_enb = num_eNodeBs;
            end
            for b_ = 1:n_enb
                [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);%UMI撒点必须在28.9范围圈内，且不能再微站的3m范围内 
                for i = 1:length(m)
                        SNR_assignment(m(i),n(i)) = 0;                    
                end
            end
            [m,n] = find(SNR_assignment(:,:)==0);
            for i = 1:length(m)
                SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));
            end
            if SYS_config.isDouble
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                for b_ = n_enb+1:num_eNodeBs
                    [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                    for i = 1:length(m)
                        SNR_assignment2(m(i),n(i)) = 0;%扇区编号等于基站编号
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
            end
        else
            % do nothing
        end
    case {'InH','InH2'}
        % do nothing
    case 'InH3'
        ISD = SYS_config.ISD;
        SNR_assignment_tmp = SNR_assignment;% 先保存原来的SNR_assignment
        SNR_assignment(:,:) = -1;% 将SNR_assignment赋成-1，-1代表不能在此处撒点的意思
        distance_matrix = networkPathlossMap.distance_matrix;
        if SYS_config.isDouble
            n_enb = num_eNodeBs/2;
        else
            n_enb = num_eNodeBs;
        end
        
        [m,n] = size(SNR_assignment);%UMA撒点的位置不能在基站的35m范围内
        for i = 1:m
            for j = 1:n
                pos=NR_common_pixel_to_pos([i j],[networkPathlossMap.roi_y(1) networkPathlossMap.roi_x(1)],networkPathlossMap.data_res);
                if pos(2)>=0 && pos(2)<=120 && pos(1)>=0 && pos(1)<=50
                    inarea = 1;
                else
                    inarea = 0;
                end
                if inarea
                    SNR_assignment(i,j) = 0; % 能够撒点的地方赋0
                end
            end
        end
        
        [m,n] = find(SNR_assignment(:,:)==0);
        for i = 1:length(m)
            SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));% 将可以撒点的地方用原来的值替换
        end
        if SYS_config.isDouble
            SNR_assignment_tmp2 = SNR_assignment2;
            SNR_assignment2(:,:) = -1;
            [m,n] = size(SNR_assignment2);
            for i = 1:m
                for j = 1:n
                    pos=NR_common_pixel_to_pos([i j],[networkPathlossMap.roi_y(1) networkPathlossMap.roi_x(1)],networkPathlossMap.data_res);
                    if pos(2)>=0 && pos(2)<=120 && pos(1)>=50+SYS_config.InH3_d && pos(1)<=100+SYS_config.InH3_d
                        inarea = 1;
                    else
                        inarea = 0;
                    end
                    if inarea
                        SNR_assignment2(i,j) = 0; % 能够撒点的地方赋0
                    end
                end
            end
            [m,n] = find(SNR_assignment2(:,:)==0);
            for i = 1:length(m)
                SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
            end
        end
    case 'UMa_to_UMi'
        if ~SYS_config.isManhattan
            % 微小区不为曼哈顿
            UMi_r = SYS_config.UMi_r;
            SNR_assignment_tmp = SNR_assignment;
            SNR_assignment(:,:) = -1;
            distance_matrix = networkPathlossMap.distance_matrix;
            micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
            for b_ = 1:num_first_sites
                [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                for i = 1:length(m)
                    SNR_assignment(m(i),n(i)) = 0;
                end
            end
            [m,n] = find(SNR_assignment(:,:)==0);
            for i = 1:length(m)
                SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));
            end
        else
            % 微小区为曼哈顿
            distance_matrix = networkPathlossMap.distance_matrix;
            O_M = [-540,-540];
            K_M = [540,540];
            O_M_pix = NR_common_pos_to_pixel(O_M, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
            K_M_pix = NR_common_pos_to_pixel(K_M, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
            M_row = O_M_pix(2):K_M_pix(2);
            M_col = O_M_pix(1):K_M_pix(1);
            tmp = zeros(size(SNR_assignment));
            tmp(M_row,M_col) = 1;
            SNR_assignment = SNR_assignment .* tmp;
        end
        
        ISD = SYS_config.ISD;
        SNR_assignment_tmp2 = SNR_assignment2;
        SNR_assignment2(:,:) = -1;
        for b_ = num_first_sites+1:num_eNodeBs
            [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
            for i = 1:length(m)
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                    if isinhex
                        SNR_assignment2(m(i),n(i)) = 0;
                    end
                end
            end
        end
        [m,n] = find(SNR_assignment2(:,:)==0);
        for i = 1:length(m)
            SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
        end
    case {'UMa_to_InH','UMi_to_InH'}
        % 求O_InH
        O_InH = eNodeB_sites(1).O_InH;
        if ~SYS_config.isFemto
            K_InH = [O_InH(1)+120,O_InH(2)+50];
        else
            K_InH = [O_InH(1)+120,O_InH(2)+40];
        end
        O_InH_pix = NR_common_pos_to_pixel(O_InH, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
        K_InH_pix = NR_common_pos_to_pixel(K_InH, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
        InH_row = O_InH_pix(2):K_InH_pix(2);
        InH_col = O_InH_pix(1):K_InH_pix(1);
        % 改inh SINR_assignment
        tmp = zeros(size(SNR_assignment));
        tmp(InH_row,InH_col) = 1;
        SNR_assignment = SNR_assignment .* tmp;
        
        switch SYS_config.scene_type
            case 'UMa_to_InH'
                ISD = SYS_config.ISD;
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                distance_matrix = networkPathlossMap.distance_matrix;
                for b_ = num_first_sites+1:num_eNodeBs
                    [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
                    for i = 1:length(m)
                        for s_ = 1:length(eNodeB_sites(b_).sectors)
                            isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                            if isinhex
                                SNR_assignment2(m(i),n(i)) = 0;
                            end
                        end
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
            case 'UMi_to_InH'
                UMi_r = SYS_config.UMi_r;
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                distance_matrix = networkPathlossMap.distance_matrix;
                micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
                for b_ = num_first_sites+1:num_eNodeBs
                    [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                    for i = 1:length(m)
                        SNR_assignment2(m(i),n(i)) = 0;
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
        end
end

% 计算扇区大小
if SYS_config.isDouble
    n_sector1 = num_first_sectors;
    n_sector2 = num_hetnet_sectors;
else
    n_sector1 = length(eNodeB_sectors);
end
cell_sizes = zeros(1,n_sector1);
cell_centers_pixel = zeros(n_sector1,2);

for s_idx = 1:n_sector1
    cell_sizes(s_idx) = sum(SNR_assignment(:)==s_idx);% 例：SNR_assignment=(1,1,2,3,3,3),sum(SNR_assignment(:)==1)得到2，意思是基站1服务了两个像素点
    [row,col] = find(SNR_assignment==s_idx);
    cell_centers_pixel(s_idx,:) = [mean(col) mean(row)];
end
if SYS_config.isDouble% 第二个系统
    cell_sizes2 = zeros(1,n_sector2);
    cell_centers_pixel2 = zeros(n_sector2,2);
    for s_idx = n_sector1+1:length(eNodeB_sectors)
        cell_sizes2(s_idx-n_sector1) = sum(SNR_assignment2(:)==s_idx);
        [row,col] = find(SNR_assignment2==s_idx);
        cell_centers_pixel2(s_idx-n_sector1,:) = [mean(col) mean(row)];
    end
    cell_centers2 = NR_common_pixel_to_pos(cell_centers_pixel2,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);% 第二个系统的扇区中心
else
    cell_sizes2 = [];
    cell_centers2 = [];
end
% 第一个系统的扇区中心
cell_centers = NR_common_pixel_to_pos(cell_centers_pixel,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);% 记录扇区中心，在画图时标记出扇区编号

%% 保存数据
networkPathlossMap.sector_assignment = SNR_assignment;
networkPathlossMap.sector_assignment_double = SNR_assignment2;
networkPathlossMap.sector_sizes = cell_sizes;% 扇区服务的总像素点数，一个1*n_sector1维数组
networkPathlossMap.sector_sizes_double = cell_sizes2;
networkPathlossMap.sector_centers = cell_centers;
networkPathlossMap.sector_centers_double = cell_centers2;

end

% 判断某点是否在六边形内
function isinhex = isInHex( dot_pix,sector_centre,ISD, networkPathlossMap)

dot_pos = NR_common_pixel_to_pos(dot_pix,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)];
for i_ = 1:7
    tmp_hex(i_,:) = sector_centre + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % 乘[1;0]表示只取第一列
end
isinhex = inpolygon(dot_pos(1),dot_pos(2),tmp_hex(:,1),tmp_hex(:,2));

end
