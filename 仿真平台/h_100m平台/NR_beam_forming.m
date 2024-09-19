function [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeBs_sectors,networkPathlossMap,SYS_config)
% 产生波束赋型图案
% 输入参数：
% UEs：用户列表
% eNodeBs_sectors：小区实体
% networkPathlossMap：路损映射矩阵
% SYS_config：LTE参数系统
%
% 输出参数：
% beam_antenna_gain：BS端波束赋型增益
% beam_ue_gain：UE端波束赋型增益

num_markers = 10; %同步情况进度条长度设置
u_markings  = round(linspace(1,length(UEs),num_markers));
DEBUG_LEVEL = SYS_config.debug_level;

ue_len=length(UEs); %所有的UE个数
enb_len=length(eNodeBs_sectors);%两个系统的话最多有122个基站（考虑UMA wrap around）,因此最大sector索引值为366
beam_antenna_gain=zeros(ue_len,enb_len );
beam_ue_gain=zeros(ue_len,enb_len );
flag=true; %用于判断每个sector的用户数是否相同
if DEBUG_LEVEL>=1
    if SYS_config.antenna_mode==2
        fprintf('Creating BS and UE beam forming: ');
    elseif SYS_config.antenna_mode==1
        fprintf('Creating BS beam forming: ');
    end
end


if SYS_config.antenna_mode==0 %确定BS端与UE端是否采用了波束赋形
    if DEBUG_LEVEL>=1
        fprintf('No BF in BS nor UE\n');
    end
else
    %导入波束赋型数据
    if SYS_config.use_cache
        try
            filename1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain.mat'];
            filename2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain.mat'];
            tmp = load(filename1);
            beam_antenna_gain = tmp.beam_antenna_gain;
            tmp = load(filename2);
            beam_ue_gain = tmp.beam_ue_gain;
            creat_beam = false;
            if DEBUG_LEVEL>=1
                fprintf('(Use beam forming cache)\n');
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
    % 初始化变量容器
    if creat_beam
        beam_antenna_gain=zeros(ue_len,enb_len);%初始化BS端波束赋型矩阵
        beam_ue_gain=zeros(ue_len,enb_len);     %初始化UE端波束赋型矩阵
        theta_map=zeros(ue_len,enb_len);        %初始化各BS到各UE的垂直角映射矩阵
        phi_map=zeros(ue_len,enb_len);          %初始化各BS到各UE的水平方位角映射矩阵
        theta_map_ue=zeros(ue_len,enb_len);     %初始化各UE到各BS的垂直角映射矩阵
        phi_map_ue=zeros(ue_len,enb_len);       %初始化各UE到各BS的水平方位角映射矩阵
        self_theta=zeros(ue_len,enb_len);       %初始化各BS在各UE被服务时BS侧波束指向的垂直角
        self_phi=zeros(ue_len,enb_len);         %初始化各BS在各UE被服务时BS侧波束指向的水平方位角
        self_phi_ue=zeros(ue_len);              %初始化各BS在各UE被服务时UE侧波束指向的垂直角
        self_theta_ue=zeros(ue_len,enb_len);    %初始化各BS在各UE被服务时UE侧波束指向的水平方位角

        for u_=1:ue_len
            if SYS_config.isDouble==false %单系统不存在邻系统干扰集合实体。故需要赋空值
                UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id=[];
            end
            
            for enb=1:enb_len
                % 如果每个扇区的UE数不等（因运动发生了切换），则flag置false
                tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id];
                if length(tmp)~=SYS_config.UE_per_eNodeB
                    flag=false;
                end
                if ~isempty (find([[UEs(u_).attached_eNodeB.in_interf_eNodeB_sectors.eNodeB_id],[UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id],UEs(u_).attached_eNodeB.eNodeB_id]==enb)) %如果b在本系统邻基站列表或者邻系统基站列表中，或者为该UE服务基站
                    switch SYS_config.scene_type
                        case {'UMA', 'UMI'} % UMA与UMI BS与UE水平角与垂直角的计算公式相同
                            % UMA/UMI垂直角的计算（sector(enb)对UE(u_)的垂直角）
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;
                            % UMA/UMI水平角的计算（sector(enb)对UE(u_)的水平角）
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            %由于水平角计算的时候减去了天线瓣的方位角，所以需要将其范围定向到[-180，180)
                            phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                            %UE(u_)对sector(enb)的方位角
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;
                            %基站水平角涉及扇区方向，UE水平角涉及UE旋转角度，因而需要将水平角的范围定向到[-180，180)
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                        case {'InH','InH2','InH3'} %InH BS的角度计算与UMA/UMI不同，UE端的角度计算与UMA/UMI相同
                            %InH BS 垂直角的计算：90-arctan(delta(y)/根号下水平距离差与站高的平方和)
                            theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                            %InH BS 水平角的计算：arctan(delta(y)/站高)
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;                          
                            %InH UE 水平角与方向角的计算与UMA/UMI相同，水平角多加了随机旋转角
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
                            % UE水平角计算的时候减去了随机的旋转角，所以需要将其范围定向到[-180，180)
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                        case 'UMa_to_UMi'%UMA干扰UMI场景角度计算方法与UMA/UMI相同，但是由于两系统是异构的，因此计算角度使用的高度差要因系统而异
                            %BS/UE水平角的计算与高度差无关，故无需改动
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;                            
                            phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                            if enb<=networkPathlossMap.num_first_sectors %基站是第一个系统的基站
                                theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            else %基站属于第二个系统
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;
                            end
                        case{'UMi_to_InH','UMa_to_InH'} %UMA/UMI对InH的干扰中，第一个系统使用InH的角度计算方法，同样要考虑异构系统的高度差得问题
                           %UE水平角的计算与其他情况相同
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                            if enb<=networkPathlossMap.num_first_sectors %基站是第一个系统的基站
                                theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
                            else %基站属于第二个系统
                                theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;                               
                            end
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                    end
                end
            end
        end
        switch SYS_config.beam_accurancy_type %判断波束偏差的类型
            case 'Constant' %偏差类型为常量
                bs_acc_th=SYS_config.beam_accuracy_theta_bs; 
                ue_acc_th=SYS_config.beam_accuracy_theta_ue;
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs;
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue;
            case 'Unique' %偏差类型为均匀分布
                bs_acc_th=SYS_config.beam_accuracy_theta_bs * rand(size(theta_map));
                ue_acc_th=SYS_config.beam_accuracy_theta_ue * rand(size(theta_map_ue));
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs * rand(size(phi_map));
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue * rand(size(phi_map_ue));
            otherwise
                error('beam accurancy type error');
        end
        %UE端朝向矩阵，找出对应第u_个UE，enb对应的服务UE的编号，从角度Map中读取，也即确定当前服务UE下所有UE的波束方向,即每个UE指向当前服务小区时，对其他小区的角度
        for u_=1:ue_len
            ser_ue_group=mod(u_,SYS_config.UE_per_eNodeB);%服务UE对应的波束组号
            if ser_ue_group==0  %如果u_可以被LTE_config.UE_per_eNodeB整除，表明其为该服务集的最后一个小区
                ser_ue_group=SYS_config.UE_per_eNodeB;
            end
            
            for enb=1:enb_len
                if ~isempty (find([[UEs(u_).attached_eNodeB.in_interf_eNodeB_sectors.eNodeB_id],[UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id],UEs(u_).attached_eNodeB.eNodeB_id]==enb)) %如果b在本系统邻基站列表或者邻系统基站列表中，或者为该UE服务基站
                    %当id为u_的UE被服务时,UE(u_)的波束指向其服务小区
                    self_theta_ue(u_,enb)=theta_map_ue(u_,UEs(u_).attached_eNodeB.eNodeB_id);
                    self_phi_ue(u_,enb)=phi_map_ue(u_,UEs(u_).attached_eNodeB.eNodeB_id);
                    if flag %如果每个小区的用户数一样，则采用分组的形式构成slot
                        %当id为u_的UE被服务时,sector(enb)的波束指向其服务ue
                        self_theta(u_,enb)=theta_map(ser_ue_group+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                        %找出对应第u_个UE，enb对应的服务UE的编号，从角度Map中读取，也即确定当前服务UE下所有sector的波束方向,即每个SECTOR指向其集合中与当前UE组号相同的本区UE
                        self_phi(u_,enb)=phi_map(ser_ue_group+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                    else %如果小区的用户数不一样，则采用小区内随机服务用户的方式
                        tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id]; %enb小区的服务UE 的id集
                        serve_ue_id=tmp(randi(length(tmp))); %从服务集中随机选取一个UE                     
                        self_theta(u_,enb)=theta_map(serve_ue_id,enb);
                        %找出对应第u_个UE，enb对应的服务UE的编号，从角度Map中读取，也即确定当前服务UE下所有sector的波束方向,即每个SECTOR指向其集合中与当前UE组号相同的本区UE
                        self_phi(u_,enb)=phi_map(serve_ue_id,enb);
                    end
                    
                end
            end
        end
        %% 4G采用36942的模型，不考虑波束赋型。NR采用38900或38901的模型，采用波束赋型
        % NR干扰LTE
        if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%如果第一个系统是LTE
            for u_ = 1:ue_len %基站侧波束赋型的计算，LTE基站对所有UE的波束赋形增益为0，NR正常计算
                if ~isempty(find(u_==u_markings,1)) %设置进度条
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
                for enb = 1:networkPathlossMap.num_first_sectors %第一个系统基站对所有UE的基站侧波束赋型为0
                    beam_antenna_gain(u_,enb)=0;
                end
                for enb = networkPathlossMap.num_first_sectors+1:enb_len
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %如果BS采用波束赋型
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%加上3dB极化增益
                    end
                end
            end
            for enb = 1:enb_len %UE侧波束赋型的计算，lteUE对所有UE的波束赋型增益为0，NR的UE正常计算
                for u_ = 1:networkPathlossMap.num_first_UEs %第一个系统UE对所有小区的UE侧波束赋形为0
                    beam_ue_gain(u_,enb)=0;
                end
                for u_ = networkPathlossMap.num_first_UEs+1:ue_len
                    if SYS_config.antenna_mode==2 %如果UE采用波束赋型
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                end
            end
            %LTE干扰NR
        elseif  ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%如果第二个系统是LTE
            for u_ = 1:ue_len %基站侧波束赋型的计算，LTE基站对所有UE的波束赋形增益为0，NR正常计算
                if ~isempty(find(u_==u_markings,1)) %设置进度条
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
                for enb = 1:networkPathlossMap.num_first_sectors %遍历第一个小区的基站
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %如果BS采用波束赋型
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%加上3dB极化增益
                    end
                end
                for enb = networkPathlossMap.num_first_sectors+1:enb_len %第二个系统基站对所有UE的基站侧波束赋型为0
                    beam_antenna_gain(u_,enb)=0;
                end
            end
            for enb = 1:enb_len %UE侧波束赋型的计算，lteUE对所有UE的波束赋型增益为0，NR的UE正常计算
                for u_ = 1:networkPathlossMap.num_first_UEs %遍历第一个小区UE
                    if SYS_config.antenna_mode==2
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                    
                end
                for u_ = networkPathlossMap.num_first_UEs+1:ue_len %第一个系统UE对所有小区的UE侧波束赋形为0
                    beam_ue_gain(u_,enb)=0;
                end
            end
            %LTE干扰LTE
        elseif strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
            if DEBUG_LEVEL>=1
                fprintf('non NR, no beam forming\n');
            end
            %NR 干扰NR
        elseif ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
            for u_=1:ue_len
                for enb=1:enb_len
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %如果BS采用波束赋型
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%加上3dB极化增益
                    end
                    if SYS_config.antenna_mode==2 %如果UE采用波束赋型
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                end
                if ~isempty(find(u_==u_markings,1)) %设置进度条
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
            end
        end
        %保存波束赋型矩阵
        if SYS_config.isSave
            filename1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain.mat'];
            filename2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain.mat'];
            save(filename1,'beam_antenna_gain');
            save(filename2,'beam_ue_gain');
        end
    end
end
end

