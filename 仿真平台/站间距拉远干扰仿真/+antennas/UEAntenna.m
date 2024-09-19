classdef UEAntenna< antennas.antenna
    %与BSAntenna功能一样，只是初始化的一些参数有所差别，当时创建时希望将UE与BS模块完全去耦合
    %因而将天线类型也分开设计，也可以更改一些基本参数只使用与BSAntenna功能一样文件
    properties
        electrical_tilt
        mechanical_tilt
        Nv
        Nh
        N_beam
        phi_3db
        theta_3db
        Am
        SLAv
        v_spacing %Vertical radiating element spacing=0.5
        h_spacing %Horizontal radiating element spacing=0.5
    end
    
    methods   
              

        
        function obj = UEAntenna(Nv,Nh,N_beam)
            obj.antenna_type = 'NR';
            %obj.max_antenna_gain =max_gain;
            obj.Nv=Nv;
            obj.Nh=Nh;
            obj.theta_3db=90;
            obj.Am=25;
            obj.SLAv=25;
            obj.N_beam=N_beam;
            obj.v_spacing=0.5;
            obj.h_spacing=0.5;
            obj.pattern_is_3D = true;
            if obj.Nh==1
%                obj.max_antenna_gain=9; %单柱/dBi
                obj.phi_3db=90;             %degree
            else
%                obj.max_antenna_gain=7.5;
                obj.phi_3db=90;
            end
                
        end
        function print(obj)
             fprintf('NR UE antenna with beam-forming\n');
        end
        
        %% 天线增益
        function h_gain=horizontal_gain(obj,phi)
            h_gain=-min(12*(phi/obj.phi_3db).^2,obj.Am);
        end
        function v_gain=vertical_gain(obj,theta)
            v_gain=-min(12*((theta-90)/obj.theta_3db).^2,obj.SLAv);
        end
        
        function Ae=elementPattern(obj,theta,phi)
            obj.max_antenna_gain=0;
% 观察max_antenna_gain的值
            Ae=obj.max_antenna_gain-min(-(obj.horizontal_gain(phi)+obj.vertical_gain(theta)),obj.Am);
            Ae = 0;
        end
  
        %% 波束赋型增益
        function antenna_gain=gain(obj,theta,phi,down_tilt_steering,horizontal_steering) %仰角，方位角，波束下倾角，波束水平角
           [a,b]=size(theta);
           temp=zeros(a,b,obj.N_beam);
           for k=1:obj.N_beam
                for n=1:obj.Nv
                    for m=1:obj.Nh                    
                    v(:,:,k,n,m)=exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*cos(theta./180.*pi)+(m-1)*obj.h_spacing.*sin(theta./180.*pi).*sin(phi./180.*pi)));
                    w(:,:,k,n,m)=1./sqrt(obj.Nh*obj.Nv).*exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*sin(down_tilt_steering(k)./180.*pi)-(m-1).*obj.h_spacing.*cos(down_tilt_steering(k)./180.*pi).*sin(horizontal_steering(k)./180.*pi)));
                    temp(:,:,k)=temp(:,:,k)+w(:,:,k,n,m).*v(:,:,k,n,m);%temp(i)为绝对值中权重值积的累和
                    end
                end
                temp(:,:,k)=abs(temp(:,:,k)).^2;
           end
           A_extra=sum(temp,3)./obj.N_beam;
          % A_extra=sum(abs(temp).^2,3)/obj.N_beam;
%           antenna_gain=obj.elementPattern(theta,phi); %+10*log10(A_extra);
%            antenna_gain=10*log10(A_extra);
            antenna_gain=obj.elementPattern(theta,phi)+10*log10(A_extra);
        end
        
        function antenna_gain=beam_gain(obj,rho,theta,phi,down_tilt_steering,horizontal_steering) %IMD3设备，只是波束增益，不包括单元增益
           [a,b]=size(theta);
           temp=zeros(a,b,obj.N_beam.^3);
          k=1;
          p=1;
          q=1;
          l=1;
             while l<=obj.N_beam
             for n=1:obj.Nv
               for m=1:obj.Nh  
                 
                         v(:,:,k,n,m)=exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*cos(theta./180.*pi)+(m-1)*obj.h_spacing.*sin(theta./180.*pi).*sin(phi./180.*pi)));                         
                         w(:,:,k,n,m)=1./sqrt(obj.Nh.*obj.Nv).*exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*sin(down_tilt_steering(p)./180.*pi)-(m-1).*obj.h_spacing.*cos(down_tilt_steering(p)./180.*pi).*sin(horizontal_steering(p)./180.*pi)))./sqrt(obj.Nh*obj.Nv).*exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*sin(down_tilt_steering(q)./180.*pi)-(m-1)*obj.h_spacing.*cos(down_tilt_steering(q)./180.*pi).*sin(horizontal_steering(q)./180.*pi))).*sqrt(obj.Nh*obj.Nv)./exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*sin(down_tilt_steering(l)./180.*pi)-(m-1)*obj.h_spacing.*cos(down_tilt_steering(l)./180.*pi).*sin(horizontal_steering(l)./180.*pi)));                         
                         temp(:,:,k)=temp(:,:,k)+w(:,:,k,n,m).*v(:,:,k,n,m);%temp(i)为绝对值中权重值积的累和
                         
               end
             end
             temp(:,:,k)=1+rho.*(abs(temp(:,:,k)).^2-1);
                          k=k+1;
                         if p==obj.N_beam
                             q=q+1;
                             p=1;
                         else
                             p=p+1;
                         end
                         if q>obj.N_beam
                             q=1;
                             l=l+1;
                         end   
             end
             A_extra=sum(temp,3)/obj.N_beam.^3;
             antenna_gain=10*log10(A_extra);
%             antenna_gain=obj.elementPattern(theta,phi);
%             antenna_gain=obj.elementPattern(theta,phi)+10*log10(A_extra);
        end
           
       function antenna_gain=gain2(obj,theta,phi,down_tilt_steering,horizontal_steering)
           obj.max_antenna_gain=12;
           obj.Am=25;
           obj.SLAv=20;
           obj.phi_3db=70;
           obj.theta_3db=10;
           Ah=-min(12*((phi-horizontal_steering)/obj.phi_3db).^2,obj.Am);
           Av=-min(12*((theta-down_tilt_steering)/obj.theta_3db).^2,obj.Am);
           antenna_gain=-min(-(Ah+Av),obj.Am);
       end
            
        function antenna_gain=gain3(obj,theta,phi,down_tilt_steering,horizontal_steering)
           [a,b]=size(theta);
           temp=zeros(a,b,obj.Nv);
            for m=1:obj.Nv
            w(:,:,m)=1/sqrt(obj.Nv).*exp(1i*2*pi*(m-1).*obj.v_spacing.*cos(down_tilt_steering./180.*pi));
            v(:,:,m)=exp(-2i*pi.*(m-1).*obj.v_spacing.*cos(theta./180.*pi));
            temp(:,:,m)=w(:,:,m).*v(:,:,m);
            end
           A_extra=abs(sum(temp,3)).^2;
            antenna_gain=10*log10(A_extra);
        end
            
        function minmaxgain = min_max_gain(obj)
            minmaxgain(1) = [0 0];
        end
    end
end
        
                
    
                    