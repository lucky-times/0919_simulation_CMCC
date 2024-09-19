classdef BSAntenna < antennas.antenna
    %���㵥Ԫ����ͼ���벨������ͼ��
    properties
        electrical_tilt % ���������
        mechanical_tilt % ��е�����
        Nv % ��ֱ�������ģ
        Nh % ˮƽ�������ģ
        N_beam % ������
        phi_3db % 3dBˮƽ����
        theta_3db % 3dB��ֱ����
        tilt %���������
        Am
        SLAv
        v_spacing %Vertical radiating element spacing=0.5
        h_spacing %Horizontal radiating element spacing=0.5
    end

    methods


        %% ���캯��
        function obj = BSAntenna(Nv,Nh,N_beam, tilt)
            obj.antenna_type = 'NR';
            %obj.max_antenna_gain =max_gain;
            obj.Nv=Nv;
            obj.Nh=Nh;
            obj.phi_3db = 15;
            obj.theta_3db = 6;
            obj.N_beam=N_beam;
            obj.v_spacing=0.5;
            obj.h_spacing=0.5;
            obj.pattern_is_3D = true;
            obj.max_antenna_gain = 20;
            obj.Am = 30;
            obj.SLAv = 30;
        end
        function print(obj)
            fprintf('NR BS antenna with beam-forming\n');
        end
%         ��Ϊ��������ģ�ͣ���ֱ����4������ˮƽ����������������������45�� 3dB���������浥������65�� 3dB����
%         % ˮƽ�������ߵ�Ԫ����
%         function h_gain=horizontal_gain(obj,phi)
%             if phi > 0
%                 h_gain=-min(12*((phi-22.5)/45).^2,obj.Am);
%             else 
%                 h_gain=-min(12*((phi-(-22.5))/45).^2,obj.Am);
%             end
%         end
%         function h_gain=horizontal_gain_65(obj,phi)
%             h_gain=-min(12*(phi/65).^2,obj.Am);
%         end
%         % ��ֱ�������ߵ�Ԫ����
%         function v_gain=vertical_gain(obj,theta)
%             v_gain=-min(12*((theta-93)/6).^2,obj.SLAv);
%         end
%         ������4���������ܵ����ߵ�Ԫͼ��%
%         function Ae=elementPattern(obj,theta,phi, tilt)
%             ˮƽ������������
%             phi2 = phi;
%             ת����ֱ�Ƕ�
%             theta_pre = 93 - theta;
% 
%             ���ƴ�ֱ����4����������-3��Ϊ��㣬��24��
%             theta = abs(theta);
%             v = ones(size(theta));
%             if abs(theta_pre)>3
%                 v = ceil((theta_pre-3)/6);
%                 v(v>3) = 3;
%                 theta2 = theta+v*6;
%             else
%                 theta2 = theta;
%             end
%             if v < 3
%                 Ae = obj.max_antenna_gain-min(-(obj.horizontal_gain(phi2)+obj.vertical_gain(theta2)),obj.Am);
%             else
%                 Ae = obj.max_antenna_gain-min(-(obj.horizontal_gain_65(phi2)+obj.vertical_gain(theta2)),obj.Am);
%             end
%         end

        %������������ģ�ͣ���ֱ����4������ˮƽ���򣬰˲�����15�� 3dB������ֱ����6�� 3dB����
        %% ˮƽ�������ߵ�Ԫ����
        function h_gain=horizontal_gain(obj,phi)
            if phi > 0
                h_gain=-min(12*((phi-7.5)/15).^2,obj.Am);
            else 
                h_gain=-min(12*((phi-(-7.5))/15).^2,obj.Am);
            end
        end

        %% ��ֱ�������ߵ�Ԫ����
        function v_gain=vertical_gain(obj,theta)
            v_gain=-min(12*((theta-(93))/6).^2,obj.SLAv);
        end
        function Ae=elementPattern(obj,theta,phi, tilt)
            %���ƴ�ֱ����4����������-6��Ϊ��㣬����3dB����24��
            alpha1 = 90 - theta;
            if alpha1>0
                v = ceil((alpha1)/6);
                if v > 3
                    v = 3;
                end
                theta2 = theta + v*6;               
            else
                theta2 = theta;
            end

            %����ˮƽ����8����������0Ϊ��㣬����3dB����120��
            if abs(phi)>7.5
                h = ceil((abs(phi)-7.5)/15);
                if h > 3
                    h = 3;
                end
                phi2 = phi - sign(phi).*h*15;
            else
                phi2 = phi;
            end
            Ae=obj.max_antenna_gain-min(-(obj.horizontal_gain(phi2)+obj.vertical_gain(theta2)),obj.Am);
        end
        %% ֻ��һ�����������ߵ�Ԫ����
        %             function Ae=elementPattern(obj,theta,phi,tilt)
        %             theta = theta + 90;
        %             theta2 = zeros(size(theta));
        %             phi2 = zeros(size(phi));
        %
        %             for i = 1:size(theta, 1)
        %                 for j = 1:size(theta, 2)
        %                     [phi2(i,j), theta2(i,j)] = utils.miscUtils.global2local(phi(i,j), theta(i,j), 0, tilt, 0);
        %                 end
        %             end
        %             theta2 = theta2 - 90;
        %
        %             Ae=obj.max_antenna_gain-min(-(obj.horizontal_gain(phi2)+obj.vertical_gain(theta2)),obj.Am);
        %         end

        %% ��������������棬������Ϊ�˶ನ������չ��ʹ����37.842��IBD3ģ�ͣ���beam_gain������
        function antenna_gain=gain(obj,theta,phi,down_tilt_steering,horizontal_steering) %���ǣ���λ�ǣ���������ǣ�����ˮƽ���
            [a,b]=size(theta);
            temp=zeros(a,b,obj.N_beam);
            for k=1:obj.N_beam
                for n=1:obj.Nv
                    for m=1:obj.Nh
                        v(:,:,k,n,m)=exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*cos(theta./180.*pi)+(m-1)*obj.h_spacing.*sin(theta./180.*pi).*sin(phi./180.*pi)));
                        w(:,:,k,n,m)=1./sqrt(obj.Nh*obj.Nv).*exp(1i.*2.*pi.*((n-1).*obj.v_spacing.*sin(down_tilt_steering(k)./180.*pi)-(m-1).*obj.h_spacing.*cos(down_tilt_steering(k)./180.*pi).*sin(horizontal_steering(k)./180.*pi)));
                        temp(:,:,k)=temp(:,:,k)+w(:,:,k,n,m).*v(:,:,k,n,m);%temp(i)Ϊ����ֵ��Ȩ��ֵ�����ۺ�
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
        %% �Ľ���IBD3��������ģ�ͣ��ڵ����������У���rho=1ʱ������ʽ������gain������ͬ
        function antenna_gain=beam_gain(obj,rho,theta,phi,down_tilt_steering,horizontal_steering) %IMD3��Ʒ
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
                        temp(:,:,k)=temp(:,:,k)+w(:,:,k,n,m).*v(:,:,k,n,m);%temp(i)Ϊ����ֵ��Ȩ��ֵ�����ۺ�

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
            antenna_gain=10*log10(A_extra);%ֻ���ز����������棬�����������뵥Ԫ����ֿ�����
            %             antenna_gain=obj.elementPattern(theta,phi);
            %              antenna_gain=obj.elementPattern(theta,phi)+10*log10(A_extra);
        end


        function minmaxgain = min_max_gain(obj)
            minmaxgain(1) = [0 0];
        end
    end
end



