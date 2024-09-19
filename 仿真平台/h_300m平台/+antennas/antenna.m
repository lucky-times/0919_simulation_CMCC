classdef antenna < handle
 % 天线模块


    properties
        antenna_type % 天线类型
        max_antenna_gain; % 最大天线增益
        pattern_is_3D = false;
    end

    methods (Abstract)

        print(obj)
        antenna_gain = gain(obj,theta)
        minmaxgain = min_max_gain(obj)
    end
    
    methods (Static)
        function attach_antenna_to_eNodeB(an_eNodeB,SYS_config,site_type)
            switch an_eNodeB.antenna_type %确定天线类型
                case 'NRAntennaBeamforming' %本平台采用使用波束赋型的扇区天线
                    switch site_type % 配置天线规模，UMA/UMI => 8*16; InH => 4*8 
                        case 'indoor'
                            an_eNodeB.antenna = antennas.BSAntenna(SYS_config.antenna_element_InH(1),SYS_config.antenna_element_InH(2),1);
                            an_eNodeB.antenna.Am=25;
                            an_eNodeB.antenna.SLAv=25;
                            an_eNodeB.antenna.theta_3db=90;
                            an_eNodeB.antenna.phi_3db=90;
                        otherwise
                            an_eNodeB.antenna = antennas.BSAntenna(SYS_config.antenna_element(1),SYS_config.antenna_element(2),1);
                    end
                otherwise
                    error('This antenna is not supported');
            end
        end
        function attach_antenna_to_UE(an_UE,SYS_config)
            switch SYS_config.antenna.antenna_gain_pattern %确定天线类型
                case 'NRAntennaBeamforming' %本平台采用使用波束赋型的扇区天线
                    % 所有场景都一样
                    an_UE.antenna=antennas.UEAntenna(2,2,1);
                    an_UE.antenna.max_antenna_gain=SYS_config.UE_max_antenna_gain;
                otherwise
                    error('This antenna is not supported');
            end
        end
    end
end
