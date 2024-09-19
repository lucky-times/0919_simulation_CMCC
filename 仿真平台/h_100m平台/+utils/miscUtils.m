classdef miscUtils
    % 包含转换角度的一些小函数
    % 另外把平台原来另一个类的两个函数也放进来了
    
    properties
    end
    
    methods(Static)
        
        % 将角度转换为[0 359]
        function mod_angle = wrapTo359(the_angle)
            mod_angle = mod(the_angle, 360);
        end
        
        % 将角度转换为[0 180]
        function mod_angle = wrapTo180(the_angle)
            the_angle(the_angle-180>0) = the_angle(the_angle-180>0)-360;
            mod_angle = the_angle;
        end
        
        % 将角度转换为[-180 180]
        function angle = wrapToAll180(angle_ini)
            % 用于将角度转换为-pi~pi
            % 输入参数：
            % angle_ini：转换前的角度（度）
            %
            % 输出参数：
            % angle：转换后的角度（度）
            
            angle_ini=angle_ini+180;
            angle_ini=mod(angle_ini,360);
            angle=angle_ini-180;
        end
        
        % 调用element类中的clear函数
        function tidy_up_memory_before_closing(UEs,eNodeBs_sectors,eNodeBs)
            for u_=1:length(UEs)
                if ~isstruct(UEs(u_))
                    UEs(u_).clear;
                end
            end
            
            for bs_=1:length(eNodeBs_sectors)
                eNodeBs_sectors(bs_).clear;
            end
            
            for site_=1:length(eNodeBs)
                eNodeBs(site_).clear;
            end
        end
        
        % 在GUI中根据listbox中选择的小区选定属于这些小区的UE
        function [UE_ids cell_sum_throughput] = get_UEs_in_given_cells(cell_ids,the_UE_traces)
            N_UEs  = length(the_UE_traces);
            cell_sum_throughput = zeros(1,length(cell_ids));
            UE_ids = false(1,N_UEs);
            for c_idx=1:length(cell_ids)
                c_ = cell_ids(c_idx);
                for u_=1:N_UEs
                    if ~isempty(find(the_UE_traces(u_).attached_eNodeB==c_,1,'first'))
                        UE_ids(u_) = true;
                    end
                end
            end
        end

        %全局坐标系转换为本地坐标系
        function [phi2,theta2]=global2local(phi1,theta1,alpha,beta,gamma)%alpha表示绕z轴旋转角度，beta表示y轴旋转角度，gamma表示x轴旋转角度,phi表示水平角，theta表示下倾角
            
            R1=[cosd(alpha),sind(alpha),0;
                -sind(alpha),cosd(alpha),0;
                0,0,1];%绕z轴旋转矩阵
            R2=[cosd(beta),0,-sind(beta);
                0,1,0;
                sind(beta),0,cosd(beta)];%绕y轴旋转矩阵
            R3=[1,0,0;
                0,cosd(gamma),sind(gamma);
                0,-sind(gamma),cosd(gamma)];%绕x轴旋转矩阵
            R = R1*R2*R3;
            direc=[sind(theta1).*cosd(phi1);
                sind(theta1).*sind(phi1);
                cosd(theta1)];
            tran=inv(R)*direc;
            theta2=acosd([0,0,1]*tran);%返回角度值
            phi2=angle([1,1j,0]*tran)/pi*180;
        end

    end
end

