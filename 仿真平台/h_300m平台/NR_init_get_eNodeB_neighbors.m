function [in_interf_eNodeB_sites, in_interf_eNodeB_sectors, ad_interf_eNodeB_sites, ad_interf_eNodeB_sectors] = NR_init_get_eNodeB_neighbors(SYS_config,this_eNodeB,eNodeB_sites,num_hetnet_sites,num_first_sites)
% 求某基站的同频及邻频干扰基站
% 输入参数: 
% SYS_config：仿真参数
% this_eNodeB：当前基站
% eNodeB_sites：所有基站
% num_hetnet_sites：异构场景中，第二个系统基站数目
% num_first_sites：第一个系统基站数
% 输出参数:   
% in_interf_eNodeB_sites：同频干扰基站
% in_interf_eNodeB_sectors：同频干扰扇区
% ad_interf_eNodeB_sites：邻频干扰基站
% ad_interf_eNodeB_sectors：邻频干扰扇区

if SYS_config.isDouble
    switch  SYS_config.scene_type
        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
            % 方法和UMa非wraparound的一样（75行）
            if this_eNodeB.id>=1 && this_eNodeB.id<=(length(eNodeB_sites)-num_hetnet_sites)
                index = 1:(length(eNodeB_sites)-num_hetnet_sites);% 第一个系统全部基站
                index(find(index(:,:)==this_eNodeB.id)) = [];% 排除自己
                in_interf_eNodeB_sites = eNodeB_sites(index);% 作为同频基站
                index = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites);% 第二个系统的全部基站作为邻频基站
                ad_interf_eNodeB_sites = eNodeB_sites(index);
            else
                index = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites);
                index(find(index(:,:)==this_eNodeB.id)) = [];
                in_interf_eNodeB_sites = eNodeB_sites(index);
                index = 1:(length(eNodeB_sites)-num_hetnet_sites);
                ad_interf_eNodeB_sites = eNodeB_sites(index);
            end
            
        otherwise
            %% 找基站
            if length(eNodeB_sites)==1
                in_interf_eNodeB_sites = [];
                ad_interf_eNodeB_sites = [];
            else
                if SYS_config.isWraparound
                    if (this_eNodeB.id>=1 && this_eNodeB.id<=num_first_sites) % 第一个系统
                        index_eNodeBs = 1:num_first_sites;% 第一个系统全部基站
                        index_eNodeBs(find(index_eNodeBs == this_eNodeB.id)) = []; % 排除掉自己
                        allDistances = zeros(1,length(index_eNodeBs));
                        for b_ = 1:length(index_eNodeBs)
                            allDistances(b_) = sqrt((eNodeB_sites(index_eNodeBs(b_)).pos(1)-this_eNodeB.pos(1))^2+(eNodeB_sites(index_eNodeBs(b_)).pos(2)-this_eNodeB.pos(2))^2);
                        end
                        [~,index] = sort(allDistances);% 假如index(1)=1指index_eNodeBs(1)距离本基站最近，index(n)=5指index_eNodeBs(5)距离本基站排第n
                        lowest_idxs = index(1:18); % 在第一个系统中找距离自己最近的18个基站
                        lowest_idxs = sort(lowest_idxs);
                        index1 = index_eNodeBs(lowest_idxs);
                        in_interf_eNodeB_sites = eNodeB_sites(index1);
                        
                        index2 = [index1 this_eNodeB.id]+num_first_sites;
                        index2 = sort(index2);
                        ad_interf_eNodeB_sites = eNodeB_sites(index2);
                    else
                        % 第二个系统基站
                        index_eNodeBs = num_first_sites+1:length(eNodeB_sites);
                        index_eNodeBs(find(index_eNodeBs == this_eNodeB.id)) = [];
                        allDistances = zeros(1,length(index_eNodeBs));
                        for b_ = 1:length(index_eNodeBs)
                            allDistances(b_) = sqrt((eNodeB_sites(index_eNodeBs(b_)).pos(1)-this_eNodeB.pos(1))^2+(eNodeB_sites(index_eNodeBs(b_)).pos(2)-this_eNodeB.pos(2))^2);
                        end
                        [~,index] = sort(allDistances);
                        lowest_idxs = index(1:18);
                        lowest_idxs = sort(lowest_idxs);
                        index1 = index_eNodeBs(lowest_idxs);
                        in_interf_eNodeB_sites = eNodeB_sites(index1);
                        
                        index2 = [index1 this_eNodeB.id]-num_first_sites;
                        index2 = sort(index2);
                        ad_interf_eNodeB_sites = eNodeB_sites(index2);
                    end
                else
                    if this_eNodeB.id>length(eNodeB_sites)/2 % 第二个系统的基站
                        index = length(eNodeB_sites)/2+1:length(eNodeB_sites);% 第二个系统的全部基站
                        index(find(index(:,:)==this_eNodeB.id)) = [];% 然后排除自己
                        in_interf_eNodeB_sites = eNodeB_sites(index);% 最终作为同频基站
                        index = 1:length(eNodeB_sites)/2;
                        ad_interf_eNodeB_sites = eNodeB_sites(index);% 第一系统的全部基站作为邻频基站
                    else % 第一个系统的基站
                        index = 1:length(eNodeB_sites)/2; % 第一个系统的全部基站
                        index(find(index(:,:)==this_eNodeB.id)) = []; % 然后排除自己
                        in_interf_eNodeB_sites = eNodeB_sites(index); % 最终作为同频基站
                        index = length(eNodeB_sites)/2+1:length(eNodeB_sites);
                        ad_interf_eNodeB_sites = eNodeB_sites(index); % 第二系统的全部基站作为邻频基站
                    end
                end
            end
    end
    
    this_eNodeB.in_interf_eNodeB_sites = in_interf_eNodeB_sites;
    this_eNodeB.ad_interf_eNodeB_sites = ad_interf_eNodeB_sites;
%             if this_eNodeB.id==62
%             fprintf('本系统干扰基站\n');
%             this_eNodeB.in_interf_eNodeB_sites.id
%             fprintf('*******\n');
%             fprintf('干扰系统干扰基站\n');
%             this_eNodeB.ad_interf_eNodeB_sites.id
%             fprintf('*******\n');
%             end
    
    %% 找扇区
    if ~isempty(in_interf_eNodeB_sites)
        for s_=1:length(this_eNodeB.sectors)
            in_interf_eNodeB_sectors = this_eNodeB.sectors([1:(s_-1) (s_+1):end]); % 在同一个基站下的同频干扰扇区
            
            for b_=1:length(in_interf_eNodeB_sites)
                in_interf_eNodeB_sectors = [in_interf_eNodeB_sectors in_interf_eNodeB_sites(b_).sectors]; % 不同基站下的同频干扰扇区
            end
            
            ad_interf_eNodeB_sectors=[];
            for b_=1:length(ad_interf_eNodeB_sites)
                ad_interf_eNodeB_sectors =[ad_interf_eNodeB_sectors ad_interf_eNodeB_sites(b_).sectors];
            end
            
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = in_interf_eNodeB_sectors;
            this_eNodeB.sectors(s_).ad_interf_eNodeB_sectors = ad_interf_eNodeB_sectors;
            %                             this_eNodeB.id
            %                             this_eNodeB.sectors(s_).eNodeB_id
            %                             fprintf('本系统干扰扇区\n');
            %                             this_eNodeB.sectors(s_).neighbors_eNodeB.eNodeB_id
            %                             fprintf('*******\n');
            %                             fprintf('干扰系统干扰扇区\n');
            %                             this_eNodeB.sectors(s_).interf_eNodeB.eNodeB_id
            %                             fprintf('*******\n');
        end
    else
        for s_=1:length(this_eNodeB.sectors)
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = [];
            this_eNodeB.sectors(s_).ad_interf_eNodeB_sectors = [];
        end
    end
    %一个系统
else
    %% 找基站
    if length(eNodeB_sites)==1
        in_interf_eNodeB_sites = [];
    else
        if SYS_config.isWraparound
            allSites = reshape([eNodeB_sites.pos],2,[])';
            allDistances = sqrt(sum((repmat(this_eNodeB.pos,[size(allSites,1) 1])-allSites).^2,2));
            [allDistances_sorted indexes] = sortrows(allDistances);
            to_take = 19;
            lowest_idxs = indexes(1:min(length(indexes),to_take));% 取离自己最近的19个基站，包括自己
            lowest_idxs = sortrows(lowest_idxs);
            lowest_idxs(find(lowest_idxs == this_eNodeB.id)) = [];
            in_interf_eNodeB_sites = eNodeB_sites(lowest_idxs);
        else
            index = 1:length(eNodeB_sites);
            index(find(index(:,:)==this_eNodeB.id)) = [];
            in_interf_eNodeB_sites = eNodeB_sites(index);
        end
    end
    
    this_eNodeB.in_interf_eNodeB_sites = in_interf_eNodeB_sites;
    %         this_eNodeB.id
    %         fprintf('!!!!!!!!\n');
    %         this_eNodeB.neighbors_eNodeB.id
    %         fprintf('*******\n');
    
    %% 找扇区
    if ~isempty(in_interf_eNodeB_sites)
        for s_=1:length(this_eNodeB.sectors)
            in_interf_eNodeB_sectors = this_eNodeB.sectors([1:(s_-1) (s_+1):end]); % 在同一个基站内的同频干扰扇区
            
            for b_=1:length(in_interf_eNodeB_sites)
                in_interf_eNodeB_sectors = [in_interf_eNodeB_sectors in_interf_eNodeB_sites(b_).sectors];% 不同基站下的同频干扰扇区
            end
            
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = in_interf_eNodeB_sectors;
        end
    else
        for s_=1:length(this_eNodeB.sectors)
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = [];
        end
    end
end
end

