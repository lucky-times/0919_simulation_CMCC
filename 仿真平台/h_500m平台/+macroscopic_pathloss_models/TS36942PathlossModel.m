classdef TS36942PathlossModel < macroscopic_pathloss_models.generalPathlossModel
     % 通过38900计算路径损耗（传播损耗）
    properties
        frequency    % 频率 MHz
        environment  % RMa和UMa选择

        Dhb          % UMa中天线高度
        Hb           % RMa中天线高度

        pathloss_function_handle % 用于调用不同场景的函数，UMa = 1 ，RMa = 2

    end

    methods
        % Class constructor
        function obj = TS36942PathlossModel(frequency,environment)
            obj.frequency = frequency;
            obj.environment = environment;
            obj.name = 'TS 36.942';
            switch environment
                case 'urban_macro'
                    obj.name = [obj.name ' urban area'];
                    obj.pathloss_function_handle = 1;
                    obj.Dhb = 25;
                case 'rural_macro'
                    obj.name = [obj.name ' rural area'];
                    obj.pathloss_function_handle = 2;
                    obj.Hb  = 45;
                otherwise
                    error(['"' environment '"" environment in TS36.942 is not valid']);
            end
        end
        
        % Returns the macroscopic pathloss in dB. Note: distance in METERS
        function [pathloss_in_db,comp] = pathloss(obj,d_2d,h_ut,d_2d_in)
            
            switch obj.pathloss_function_handle
                case 1
                    d_3d=sqrt(d_2d.^2+(obj.Dhb-h_ut).^2);
                    d_3d_in = d_3d.*d_2d_in./d_2d;
                    d_3d_out = d_3d-d_3d_in;
                    pathloss_in_db = obj.pathloss_urban(d_3d_out);
                case 2
                    d_3d=sqrt(d_2d.^2+(obj.Hb-h_ut).^2);
                    d_3d_in = d_3d.*d_2d_in./d_2d;
                    d_3d_out = d_3d-d_3d_in;
                    pathloss_in_db = obj.pathloss_rural(d_3d_out);
            end
            comp = zeros(size(d_2d));% 都只适用于NLOS
        end
        
        % Urban area pathloss
        function pl = pathloss_urban(obj,distance)
            %           distance ... actual distance in m
            % output:   pl       ... NLOS pathloss in dB
            
            distance = distance/1e3;          % Calculations are done in Km
            frequency = obj.frequency/1e9; % Calculations are done in freq in GHz
            pl = 28 + 22*log10(distance) + 20*log10(frequency);
%             pl = 40*(1-4e-3*obj.Dhb)*log10(distance)-18*log10(obj.Dhb)+21*log10(frequency)+80;
        end
        
        % Rural area pathloss
        function pl = pathloss_rural(obj,distance)
            %           distance ... actual distance in m
            % output:   pl       ... NLOS pathloss in dB
            
            distance = distance/1000;          % Calculations are done in Km
            frequency = obj.frequency/1000000; % Calculations are done in freq in MHz

            pl = 69.55+26.16*log10(frequency)-13.82*log10(obj.Hb)+(44.9-6.55*log10(obj.Hb))*log10(distance)-4.78*(log10(frequency)^2)+18.33*log10(frequency)-40.94;
        end
    end
end