classdef downlinkChannelModel < handle
% 指定UE的下行信道模型——每个UE须有一个特定的信道模型，包含一些链路信息。
% 之所以叫downlinkChannelModel，是因为路损、阴影衰落、穿透损耗等损耗信息对于上下行链路来说是相等的
% 该类的对象作为UE的一个属性保存在UE中，即一个UE就有一条链路信息

   properties

       macroscopic_pathloss_model_is_set = false;
       macroscopic_pathloss_model % 路损模型
       shadow_fading_model_is_set = false; 
       shadow_fading_model % 阴影衰落模型 
       attached_UE % 该信道归属的UE
   end

   methods
       %% 构造函数
       function obj = downlinkChannelModel(aUE)
           obj.attached_UE = aUE;
       end
       
       %% 下面一系列函数的作用是调用macroscopicPathlossMap或shadowFadingMapClaussen的函数获取相关数据
       % 大体方法是：通过attached_UE调用UE中的id，pos，所属的小区，进而作为参数调用macroscopicPathlossMap或shadowFadingMapClaussen
       % 中的参数。有些求干扰径的损耗的函数需要调用一系列的UE（如interferingEnodeBids是一个1*N维的矩阵），所以与服务径的分开写
       
       %% 得到UE的高度
         function height = UE_height(obj)
           pos = obj.attached_UE.pos;
           height = obj.macroscopic_pathloss_model.get_UE_height(pos);
         end
         %% 得到UE的旋转角
         function orientation = UE_orientation(obj)
           pos = obj.attached_UE.pos;
           orientation = obj.macroscopic_pathloss_model.get_UE_orientation(pos);
         end
        %% 得到UE到其服务小区的穿透损耗 
       function penetration_loss = macroscopic_penetration_loss(obj)
           eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
           pos = obj.attached_UE.pos;
           penetration_loss = obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(pos,eNodeB_id);
       end
       %% 得到UE到某个干扰小区的穿透损耗
       function penetration_loss= interfering_penetration_loss(obj,interferingEnodeBids)
           pos = obj.attached_UE.pos;
           penetration_loss = reshape(obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(pos,interferingEnodeBids),[],1);
       end
       %% 得到上行干扰UE到该UE服务小区的穿透损耗
       function penetration_loss= ul_interfering_penetration_loss(obj,this_sector_id,interfering_ue_id)
           UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
           interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
           for u_ = 1:size(interfering_ue_pos,1);
               penetration_loss(u_) = obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),this_sector_id);
           end
       end
         %% 得到该UE到服务小区的传播损耗
       function pathloss_900 = macroscopic_pathloss_900(obj)
           eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
           pos = obj.attached_UE.pos;
           pathloss_900 = obj.macroscopic_pathloss_model.get_pathloss_900_eNodeB(pos,eNodeB_id);
       end
       
       %% 得到指定干扰小区到该UE的传播损耗
       function pathloss_900=  interfering_macroscopic_pathloss_900(obj,interferingEnodeBids)
           pos = obj.attached_UE.pos;
           pathloss_900 = reshape(obj.macroscopic_pathloss_model.get_pathloss_900_eNodeB(pos,interferingEnodeBids),[],1);
      end
       %% 得到服务小区到该UE的传播损耗
       function [pathloss is_set] = macroscopic_pathloss(obj)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else

               eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
               pos = obj.attached_UE.pos;
               pathloss = obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,eNodeB_id);
               is_set   = true;
           end
       end
       
       %% 得到指定干扰小区到该UE的传播损耗
       function [pathloss is_set] = interfering_macroscopic_pathloss(obj,interferingEnodeBids)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = zeros([length(interferingEnodeBids) 1]); % Column vector
               is_set   = false;
               return
           else
               pos = obj.attached_UE.pos;
               pathloss = reshape(obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,interferingEnodeBids),[],1);
               is_set   = true;
           end
       end
       %% 得到上行干扰UE到该服务基站的路损
       function [pathloss is_set] = ul_interfering_macroscopic_pathloss(obj,this_sector_id,interfering_ue_id)
           if ~obj.macroscopic_pathloss_model_is_set
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               pathloss = zeros(1,size(interfering_ue_pos,1));
               is_set   = false;
               return
           else
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               for u_ = 1:size(interfering_ue_pos,1);
                   pathloss(u_) = obj.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),this_sector_id);
               end
               is_set   = true;
           end
       end
       
       %% 得到该UE的阴影衰落
       function [pathloss is_set] = shadow_fading_pathloss(obj)
           if ~obj.shadow_fading_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else
               % 得到eNodeB id 和pos
               eNodeB_id = obj.attached_UE.attached_site.id;
               pos       = obj.attached_UE.pos;
               pathloss  = obj.shadow_fading_model.get_pathloss(pos,eNodeB_id);
               is_set    = true;
           end
       end
       
       %% 得到指定干扰基站到该UE的阴影衰落
       function [pathloss is_set] = interfering_shadow_fading_pathloss(obj,interferingSiteIds)
           if ~obj.shadow_fading_model_is_set
               pathloss = zeros([length(interferingSiteIds) 1]); % Column vector
               is_set   = false;
               return
           else
               pos      = obj.attached_UE.pos;
               pathloss = reshape(obj.shadow_fading_model.get_pathloss(pos,interferingSiteIds),[],1);
               is_set   = true;
           end
       end
       %% 得到干扰UE到该服务基站的阴影衰落
       function [pathloss is_set] = ul_interfering_shadow_fading_pathloss(obj,this_enb_id,interfering_ue_id)
           if ~obj.shadow_fading_model_is_set
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               pathloss = zeros(1,size(interfering_ue_pos,1));
               is_set   = false;
               return
           else
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               for u_ = 1:size(interfering_ue_pos,1);
                   pathloss(u_) = obj.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),this_enb_id);
               end
               is_set   = true;
           end
       end
       %% 得到RB分配情况表
       function the_RB_grid = RB_grid(obj)
           the_RB_grid = obj.attached_UE.attached_eNodeB.RB_grid;
       end
       %% 设置宏路损模型
       function set_macroscopic_pathloss_model(obj,macroscopic_pathloss_model)
           obj.macroscopic_pathloss_model        = macroscopic_pathloss_model;
           obj.macroscopic_pathloss_model_is_set = true;
       end
       %% 设置阴影衰落模型
       function set_shadow_fading_model(obj,shadow_fading_model)
           obj.shadow_fading_model        = shadow_fading_model;
           obj.shadow_fading_model_is_set = true;
       end
   end
end 
