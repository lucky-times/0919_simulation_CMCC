function [ue_or_bs_gain_asy,ue_gain_asy,bs_or_ue_gain_asy,bs_gain_asy] = NR_beam_forming_asy(UEs,eNodeBs_sectors,networkPathlossMap,SYS_config)
% 计算不同步场景下的波束赋型矩阵（一定需要双系统）
% 产生波束赋型图案
% 输入参数：
% UEs：用户列表
% eNodeBs_sectors：小区实体
% networkPathlossMap：路损映射矩阵
% SYS_config：LTE参数系统
%
% 输出参数：
% ue_or_bs_gain_asy：上行干扰下行中，发射端波束赋型矩阵——本系统系统内BS端波束赋型增益，邻系统中对应为干扰UE波束赋型增益
% ue_gain_asy：上行干扰下行中，接收端波束赋型矩阵——UE端波束赋型增益
% bs_or_ue_gain_asy：下行干扰上行中，发射端波束赋型矩阵——本系统内UE波束赋型增益，邻系统中对应为干扰BS波束赋型增益
% bs_gain_asy： 下行干扰上行中，接收端波束赋型矩阵——BS端波束赋型增益
flag=true; %用于标识两系统每BS的UE数是否相同
DEBUG_LEVEL = SYS_config.debug_level;
num_markers = 10; %进度条长度
u_markings  = round(linspace(1,length(UEs),num_markers));
ue_len=length(UEs);
enb_len=length(eNodeBs_sectors);%两个系统的话最多有122个基站（考虑UMA wrap around）,因此最大sector索引值为366
ue_or_bs_gain_asy=zeros(ue_len,enb_len );
ue_gain_asy=zeros(ue_len,enb_len );
bs_or_ue_gain_asy=zeros(ue_len,enb_len );
bs_gain_asy=zeros(ue_len,enb_len );

if DEBUG_LEVEL>=1
    if SYS_config.antenna_mode==2
        fprintf('Creating BS and UE beam forming for asynchronization:');
    elseif SYS_config.antenna_mode==1
        fprintf('Creating BS beam forming for asynchronization:');
    end
end


if SYS_config.antenna_mode==0
    if DEBUG_LEVEL>=1
        fprintf('No BF in BS nor UE\n')
    end
else
    %导入波束赋型数据
    if SYS_config.use_cache
        file_name1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_up_down.mat'];
        file_name2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_up_down.mat'];
        file_name3 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_down_up.mat'];
        file_name4 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_down_up.mat'];
        try
            tmp = load(file_name1);
            ue_or_bs_gain_asy = tmp.ue_or_bs_gain_asy;
            tmp = load(file_name2);
            ue_gain_asy = tmp.ue_gain_asy;
            tmp = load(file_name3);
            bs_gain_asy = tmp.bs_gain_asy;
            tmp = load(file_name4);
            bs_or_ue_gain_asy = tmp.bs_or_ue_gain_asy;
            creat_beam = false;
            if DEBUG_LEVEL>=1
                fprintf('(use beam forming cache)\n');
            end
        catch
            if DEBUG_LEVEL>=1
                fprintf('Can not find %s,creating it:',SYS_config.beam_cache_file);
                fprintf('\n               ');
            end
            creat_beam = true;
        end
    else
        creat_beam = true;
    end
    %% 初始化容器
    if creat_beam
        ue_or_bs_gain_asy=zeros(ue_len,enb_len);     %初始化上行干扰下行时发射端波束赋形矩阵（本系统BS发出信号，邻系统UE发出信号）
        ue_gain_asy=zeros(ue_len,enb_len);           %初始化上行干扰下行时接收端波束赋型矩阵（本系统与邻系统UE端接收信号）
        bs_or_ue_gain_asy=zeros(ue_len,enb_len);     %初始化下行干扰上行时发射端波束赋型矩阵（本系统UE发出信号，邻系统BS发出信号）
        bs_gain_asy=zeros(ue_len,enb_len);           %初始化下行干扰上行时接收端波束赋形矩阵（本系统与邻系统BS端接收信号）
        theta_map_ue=zeros(ue_len,enb_len);          %初始化各UE对其所在系统各基站的垂直角矩阵（系统内干扰形式不变）
        phi_map_ue=zeros(ue_len,enb_len);            %初始化各UE对其所在系统各基站的水平角矩阵
        theta_map=zeros(ue_len,enb_len);             %初始化各BS对其所在系统各UE的垂直角矩阵
        phi_map=zeros(ue_len,enb_len);               %初始化各BS对其所在系统各UE的水平角矩阵
        self_theta=zeros(ue_len,enb_len);
        self_theta_ue=zeros(ue_len,enb_len);
        self_phi=zeros(ue_len,enb_len);
        self_phi_ue=zeros(ue_len);
        beam_orientation_phi_ue2bs=zeros(ue_len);    %初始化UE对其服务BS的波束水平方向角
        beam_orientation_phi_bs2ue=zeros(ue_len);    %初始化服务BS对UE的波束水平方向角
        s=zeros(ue_len);                             %UE服务基站的id
        beam_orientation_theta_ue2bs=zeros(ue_len);  %初始化UE对其服务BS的波束方向垂直角
        beam_orientation_theta_bs2ue=zeros(ue_len);  %初始化服务BS对UE的波束垂直角
        ser_group_num=zeros(ue_len);                 %快照组号
        h1=zeros(ue_len);                            %初始化vic系统基站的高度
        h2=zeros(ue_len);                            %初始化agg系统基站的高度
        phi_up_down_agg=zeros(ue_len,enb_len);       %初始化agg系统UE对vic系统UE的水平角
        phi_up_down_vic=zeros(ue_len,enb_len);       %初始化vic系统UE对agg系统UE的水平角
        theta_up_down_agg=zeros(ue_len,enb_len);     %初始化agg系统UE对vic系统UE的垂直角
        theta_up_down_vic=zeros(ue_len,enb_len);     %初始化vic系统UE对agg系统UE的垂直角
        phi_down_up_agg=zeros(ue_len,enb_len);       %初始化agg系统BS对vic系统BS的水平角
        phi_down_up_vic=zeros(ue_len,enb_len);       %初始化vic系统BS对agg系统BS的水平角
        theta_down_up_agg=zeros(ue_len,enb_len);     %初始化agg系统BS对vic系统BS的垂直角
        theta_down_up_vic=zeros(ue_len,enb_len);     %初始化vic系统BS对agg系统BS的垂直角
        %对每个UE配置UE的天线，并计算每个UE的波束方向
        for enb=1:enb_len
            % 如果每个扇区的UE数不等（因运动发生了切换），则flag置false
            tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id];
            if length(tmp)~=SYS_config.UE_per_eNodeB
                flag=false;
            end
        end
        
        for u_=1:ue_len
            s(u_)=UEs(u_).attached_eNodeB.eNodeB_id; %用来确定改UE服务的sector
            beam_orientation_phi_ue2bs(u_)=atan2(eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;%计算UEs(u_)对其服务基站的水平角
            %如果场景为InH或者为UMa_to_InH、UMi_to_InH场景下第一个系统时，采用InH的角度计算方法
            if strcmp ('InH',SYS_config.scene_type) || strcmp ('InH2',SYS_config.scene_type) ||strcmp ('InH3',SYS_config.scene_type) ||( u_<= networkPathlossMap.num_first_UEs && (strcmp('UMa_to_InH',SYS_config.scene_type) || strcmp('UMi_to_InH',SYS_config.scene_type)))
                beam_orientation_phi_bs2ue(u_)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
            else
                beam_orientation_phi_bs2ue(u_)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180-UEs(u_).orientation;%计算其服务基站对UE(u_)的水平角
            end
            % 将服务水平角转换为(-180,180]
            beam_orientation_phi_ue2bs(u_)=utils.miscUtils.wrapToAll180(beam_orientation_phi_ue2bs(u_));
            beam_orientation_phi_bs2ue(u_)=utils.miscUtils.wrapToAll180(beam_orientation_phi_bs2ue(u_));
            switch SYS_config.scene_type
                case {'InH','InH2','InH3'}
                    beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                    beam_orientation_theta_bs2ue(u_)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180;%计算服务基站对UE(u_)的垂直角
                case {'UMa_to_InH','UMi_to_InH'}
                    if u_<= networkPathlossMap.num_first_UEs %如果是第一个系统，按照InH方法计算服务角度
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180;%计算服务基站对UE(u_)的垂直角
                    else %如果是第二个系统
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;%计算UE(u_)对其服务基站的垂直角
                    end
                case 'UMa_to_UMi'
                    if u_<= networkPathlossMap.num_first_UEs %UE来自第一个系统
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                    else %UE来自第二个系统
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;%计算UE(u_)对其服务基站的垂直角
                    end
                otherwise
                    beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%计算UE(u_)对其服务基站的垂直角
                    beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;%计算UE(u_)对其服务基站的垂直角

            end
            
            if flag % 如果每个sector的用户数相同则按序计算UEs(u_)的服务组号
                ser_group_num(u_)=mod(u_,SYS_config.UE_per_eNodeB);
                if ser_group_num(u_)==0 %如果取余为0，服务组号设为LTE_config.UE_per_eNodeB
                    ser_group_num(u_)=SYS_config.UE_per_eNodeB;
                end
            else %如果每个sector用户数不同，则随机在每个用户的UE列表中选出一个
                tmp= [UEs(u_).attached_eNodeB.attached_UEs_vector.id];
                rand_index=randi(length(tmp));
                serve_id(u_)=tmp(rand_index);
            end
            
            
        end
        %计算UE相对于另一个UE的波束赋型
        for u_=1:ue_len
            for enb=1:enb_len
                if (u_<= networkPathlossMap.num_first_UEs &&enb>networkPathlossMap.num_first_sectors) || (u_>networkPathlossMap.num_first_UEs && enb<=networkPathlossMap.num_first_sectors)
                    if enb>networkPathlossMap.num_first_sectors %基站来自第二个系统
                       h1(u_)=SYS_config.site_height; %h1表示vic系统基站的高度
                       switch SYS_config.scene_type
                           case {'UMa_to_UMi','UMa_to_InH', 'UMi_to_InH'}
                               h2(u_)=SYS_config.site_height2; %h2表示agg系统基站高度
                           otherwise
                               h2(u_)=SYS_config.site_height;
                       end
                   else %干扰基站来自第一个系统，第二个系统是vic系统
                       h2(u_)=SYS_config.site_height;
                       switch SYS_config.scene_type
                           case {'UMa_to_UMi','UMa_to_InH', 'UMi_to_InH'}
                               h1(u_)=SYS_config.site_height2;
                           otherwise
                               h1(u_)=SYS_config.site_height;
                       end
                   end
                   
                if flag                                       
                    %计算agg系统UE对vic系统UE的水平角(上行干扰下行,从干扰UE的发射增益出发)
                    phi_up_down_agg(u_,enb)=atan2(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2),UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1))./pi*180-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).orientation;%干扰UE对UE(u_)发射赋型增益
                    %计算vic系统UE对agg系统UE的水平角(上行干扰下行,从被干扰UE的接收增益出发)
                    phi_up_down_vic(u_,enb)=atan2(UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)-UEs(u_).pos(2),UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    %计算agg系统UE对vic系统UE的垂直角
                    theta_up_down_agg(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)).^2),UEs(u_).height-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height)./pi*180;
                    %计算vic系统UE对agg系统UE的垂直角
                    theta_up_down_vic(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)).^2),UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height-UEs(u_).height)./pi*180;
                else 
                    phi_up_down_agg(u_,enb)=atan2(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2),UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1))./pi*180-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).orientation;%干扰UE对UE(u_)发射赋型增益
                    phi_up_down_vic(u_,enb)=atan2(UEs(serve_id(u_)).pos(2)-UEs(u_).pos(2),UEs(serve_id(u_)).pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    theta_up_down_agg(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2)).^2),UEs(u_).height-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height)./pi*180;
                    theta_up_down_vic(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2)).^2),UEs(serve_id).height-UEs(u_).height)./pi*180;
                end

                    %计算agg系统sector对vic系统sector的水平角（下行干扰上行,从干扰sector的发射增益出发）
                    phi_down_up_agg(u_,enb)=atan2(UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))-eNodeBs_sectors(enb).azimuth;
                    %计算vic系统sector对agg系统sector的水平角（下行干扰上行,从被干扰sector的接收增益出发）
                    phi_down_up_vic(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(1))-UEs(u_).attached_eNodeB.azimuth;
                    %角度转换为(-180,180]
                    phi_up_down_agg(u_,enb)=utils.miscUtils.wrapToAll180(phi_up_down_agg(u_,enb));
                    phi_up_down_vic(u_,enb)=utils.miscUtils.wrapToAll180(phi_up_down_vic(u_,enb));
                    phi_down_up_agg(u_,enb)=utils.miscUtils.wrapToAll180(phi_down_up_agg(u_,enb));
                    phi_down_up_vic(u_,enb)=utils.miscUtils.wrapToAll180(phi_down_up_vic(u_,enb));
    
                    %计算agg系统sector对vic系统sector的垂直角
                    theta_down_up_agg(u_,enb)=atan2(sqrt((UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),h1(u_)-h2(u_))./pi*180;
                    %计算vic系统sector对vic系统sector的垂直角
                    theta_down_up_vic(u_,enb)=atan2(sqrt((eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)).^2+(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)).^2),h2(u_)-h1(u_))./pi*180;

                else
                    %UE端考虑水平旋转角，并且垂直角与水平角与基站端的互补
                    phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),h1(u_)-UEs(u_).height)./pi*180;
                    
                     %h1记录的是victim的系统站高，而当遍历到第一个系统的UE时，第一个系统为victim系统，遍历到第二个系统的UE时，第二个系统为相应的"victim系统"
                     %这里相当于令h1为u_所在系统的站高。
                    switch SYS_config.scene_type
                        case {'UMA', 'UMI','RMa','UMa_to_UMi'}
                            % UMA/UMI垂直角的计算
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-h1(u_))./pi*180;
                            % UMA/UMI水平角的计算：arctan(delta(y)/delta(x))-sector azimuth
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                        case {'InH','InH2','InH3'}
                            %InH BS 垂直角的计算：90-arctan(delta(y)/根号下水平距离差与站高的平方和)
                            theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                            %InH BS 水平角的计算：arctan(delta(y)/站高)
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),h1(u_)-SYS_config.UE_height)./pi*180;                            
                        case {'UMa_to_InH','UMi_to_InH'}
                            if enb <= networkPathlossMap.num_first_sectors %第一个系统采用InH的方式
                                theta_map(u_,enb)=atan2(sqrt((h1(u_)- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),h1(u_)-UEs(u_).height)./pi*180;
                            else %第二个系统采用UMA形式
                                theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-h1(u_))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            end
                                
                    end
                    % UE水平角计算的时候减去了随机的旋转角，所以需要将其范围定向到[-180，180)
                    phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                    %由于水平角计算的时候减去了天线瓣的方位角，所以需要将其范围定向到[-180，180)
                    phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                end
            end
        end
        
        switch SYS_config.beam_accurancy_type %判断波束偏差的类型
            case 'Constant' %偏差类型为常量
                bs_acc_th=SYS_config.beam_accuracy_theta_bs;
                ue_acc_th=SYS_config.beam_accuracy_theta_ue;
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs;
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue;
            case 'Unique'  %偏差类型为均匀分布
                bs_acc_th=SYS_config.beam_accuracy_theta_bs * rand(size(theta_map));
                ue_acc_th=SYS_config.beam_accuracy_theta_ue * rand(size(theta_map_ue));
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs * rand(size(phi_map));
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue * rand(size(phi_map_ue));
            otherwise
                error('beam accurancy type error')
        end
        
        for u_=1:ue_len
            for enb=1:enb_len
                if (u_<= networkPathlossMap.num_first_UEs &&enb>networkPathlossMap.num_first_sectors) || (u_>networkPathlossMap.num_first_UEs && enb<=networkPathlossMap.num_first_sectors)
                    if SYS_config.antenna_mode==2
                        %（上行干扰下行）UE(u_)的接收赋型增益（相对邻系统）
                        ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_up_down_vic(u_,enb),phi_up_down_vic(u_,enb),beam_orientation_theta_ue2bs(u_)-90+ue_acc_th,beam_orientation_phi_ue2bs(u_)+ue_acc_phi);
                        %（下行干扰上行）UE(u_)归属的sector接收赋型增益（相对邻系统）
                        bs_gain_asy(u_,enb) = UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_down_up_vic(u_,enb),phi_down_up_vic(u_,enb),beam_orientation_theta_bs2ue(u_)-90+bs_acc_th,beam_orientation_phi_bs2ue(u_)+bs_acc_phi);
                        if flag
                            %（上行干扰下行）UE(u_)受到的UE（邻系统）干扰发射赋型增益
                            ue_or_bs_gain_asy(u_,enb)=UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).antenna.beam_gain(1,theta_up_down_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),phi_up_down_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),beam_orientation_theta_ue2bs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))-90+ue_acc_th,beam_orientation_phi_ue2bs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))+ue_acc_phi);
                            %（下行干扰上行）UE(u_)归属的sector受到sector（邻系统）干扰发射赋型增益
                            bs_or_ue_gain_asy(u_,enb) = eNodeBs_sectors(enb).antenna.beam_gain(1,theta_down_up_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),phi_down_up_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),beam_orientation_theta_bs2ue((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))-90+bs_acc_th,beam_orientation_phi_bs2ue((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))+bs_acc_phi);
                        else
                            ue_or_bs_gain_asy(u_,enb)=UEs(serve(u_)).antenna.beam_gain(1,theta_up_down_agg(serve_id(u_),s(u_)),phi_up_down_agg(serve_id(u_),s(u_)),beam_orientation_theta_ue2bs(serve_id(u_))-90+ue_acc_th,beam_orientation_phi_ue2bs(serve_id(u_))+ue_acc_phi);
                            bs_or_ue_gain_asy(u_,enb) = eNodeBs_sectors(enb).antenna.beam_gain(1,theta_down_up_agg(serve_id(u_),s(u_)),phi_down_up_agg(serve_id(u_),s(u_)),beam_orientation_theta_bs2ue(serve_id(u_))-90+bs_acc_th,beam_orientation_phi_bs2ue(serve_id(u_))+bs_acc_phi);
                        end
                    end
                 else
                    % 计算当前UE（u_）被服务时，enb与其服务UE所成的角
                    if ~flag %如果两系统每个基站服务的用户数不一样，则随机选择服务UE取快照
                        self_theta(u_,enb)=theta_map(serve_id(u_),enb);
                        self_phi(u_,enb)=phi_map(serve_id(u_),enb);
                    else %如果两系统每个基站服务的用户数一样，则按序选择服务UE取快照
                        self_theta(u_,enb)=theta_map(ser_group_num(u_)+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                        self_phi(u_,enb)=phi_map(ser_group_num(u_)+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                    end
                    self_theta_ue(u_,enb)=theta_map_ue(u_,s(u_));
                    self_phi_ue(u_,enb)=phi_map_ue(u_,s(u_));
                    if SYS_config.antenna_mode==2
                        ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                        bs_or_ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2
                        ue_or_bs_gain_asy(u_,enb)= UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%加上3dB极化增益
                        bs_gain_asy(u_,enb)= UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%加上3dB极化增益
                    end
                end
            end
            if ~isempty(find(u_==u_markings,1)) %设置进度条
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
        end
    end
    %NR干扰LTE
    if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%如果第一个系统是LTE
        for u_=1:ue_len
            if ~isempty(find(u_==u_markings,1)) %设置进度条
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
            for enb = 1:networkPathlossMap.num_first_sectors   %第一个系统基站对第一个系统的UE的基站侧赋型增益为0，第一个系统UE对第二个系统的UE的接收赋型增益为0
                ue_or_bs_gain_asy(u_,enb)=0;
                bs_or_ue_gain_asy(u_,enb)=0;
            end
        end
        for enb=1:enb_len
            for u_=1:networkPathlossMap.num_first_UEs %第一个系统UE对第一个系统的基站的UE侧赋型增益为0，第一个系统UE对第二个系统的发射赋型增益为0
                ue_gain_asy(u_,enb)=0;
                bs_gain_asy(u_,enb)=0;
            end
        end
        %LTE干扰NR
    elseif ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%如果第二个系统是LTE
        for u_=1:ue_len
            if ~isempty(find(u_==u_markings,1)) %设置进度条
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
            for enb = networkPathlossMap.num_first_sectors+1:enb_len   %第二个系统基站对第二个系统的UE的基站侧赋型增益为0，第二个系统UE对第一个系统的UE的接收赋型增益为0
                ue_or_bs_gain_asy(u_,enb)=0;
                bs_or_ue_gain_asy(u_,enb)=0;
            end
        end
        for enb=1:enb_len
            for u_=networkPathlossMap.num_first_UEs+1:ue_len    %第二个系统UE对第二个系统的基站的UE侧赋型增益为0，第二个系统UE对第一个系统的发射赋型增益为0
                ue_gain_asy(u_,enb)=0;
                bs_gain_asy(u_,enb)=0;
            end
        end
    elseif strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
        if DEBUG_LEVEL>=1
            fprintf('non NR, no beam forming\n');
        end
        ue_or_bs_gain_asy=zeros(ue_len,enb_len);
        ue_gain_asy=zeros(ue_len,enb_len);
    end
    
    %保存波束赋型矩阵
    if SYS_config.isSave
        file_name1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_up_down.mat'];
        file_name2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_up_down.mat'];
        file_name3 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_down_up.mat'];
        file_name4 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_down_up.mat'];
        save(file_name1,'ue_or_bs_gain_asy');
        save(file_name2,'ue_gain_asy');
        save(file_name3,'bs_gain_asy');
        save(file_name4,'bs_or_ue_gain_asy');
    end
end

end
