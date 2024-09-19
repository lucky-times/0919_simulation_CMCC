function eNodeBs = LTE_init_read_atoll_cell_data(capesso_params)
% Reads out the specified folder and searchs for the following files:
%   - *.par files 
%   - LTE_BTS_trial_cluster_20100809.txt
%   - LTE_sites_trial_cluster.txt
%   - LTE_cell_trial_cluster.txt
%
%   Text files have a different format in newer Atoll version
%   For files of old Atoll version use function
%   LTE_init_read_atoll_cell_data
% 
% (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010

%% Parameter definitions
folder                        = capesso_params.pathloss_data_folder;
cell_atoll_filename           = capesso_params.cell_atoll_filename;
site_atoll_filename           = capesso_params.site_atoll_filename;
transmitter_atoll_filename    = capesso_params.transmitter_atoll_filename;
kathrein_antenna_folder       = capesso_params.kathrein_antenna_folder;

% Workaround for using frequency defined in LTE_load_params
if capesso_params.frequency == 2.14e+9
    frequency = 2140;
else
    warning('Selected frequency for the antenna gain pattern has been substituted for convinience for 2.14 GHz');
end

%% Read out LTE_cell_trial_cluster
% The fields are
% - Name	
% - Transmitter	
% - Max Power (dBm)	
% - Frequency Band	
% - Channel Number	
% - Traffic Load (UL) (%)	
% - Traffic Load (DL) (%)	
% - LTE Equipment	
% - Comments	
% - Max Number of Users	
% - UL Noise Rise (dB)	
% - Max number of intra-technology neighbours	
% - Max number of inter-technology neighbours	
% - Active	AMS & MU-MIMO Threshold (dB)
% - Scheduler	
% - Reference Signal C/N Threshold (dB)	
% - Diversity Support (DL)	
% - TDD Frame Configuration	
% - Physical Cell ID	
% - Min Reuse Distance (m)	
% - Physical Cell ID Status	
% - Channel Allocation Status	
% - Max Traffic Load (UL) (%)
% - Max Traffic Load (DL) (%)
% - Inter-technology UL Noise Rise (dB)
% - Inter-technology DL Noise Rise (dB)
% - Diversity Support (UL)
% - MU-MIMO Capacity Gain (UL)
% - SS & PBCH EPRE Offset / RS (dB)
% - PDSCH & PDCCH EPRE Offset / RS (dB)
% - Layer (Lowest Layer = Highest Priority)
% - Interference Coordination Support
% - ICIC Ratio (DL) (%)
% - ICIC Delta Path Loss Threshold (dB)
% - Fractional Power Control Factor
% - Max UL Noise Rise (dB)
% - Max PUSCH C/(I+N) (dB)
% - ICIC UL Noise Rise (dB)
% - RS_EPRE	PBCH_POWER_OFFSET
% - PDCCH_POWER_OFFSET
% - PSS ID
% - SSS ID
% - Instantaneous Reference Signal Power (dBm)
% - Instantaneous SS & PBCH Power (dBm)
% - Average PDSCH & PDCCH Power (dBm)	


% Since the numerical fields are separated by a COMMA instead of a point, I
% will have to do the 4,7 -> 4.7 conversion by hand
fid_ = fopen (fullfile(folder,cell_atoll_filename),'r');
% M = textscan(fid_,'%s %s %n %s %n %n %n %n %n %s %n %n %s %n %n %n %n %n %s %s %n %n %n %n %s %s %n %n %n %s %s %n %n %s %s %n %n %n %n %n %n','Delimiter','\t','Headerlines',1);
cell_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
fclose(fid_);

%% Read LTE_sites_trial_cluster

% The fields are:
% - Name
% - X
% - Y
% - Altitude (m)
% - Pylon Height (m)
% - Support Type
% - SITENAME
% - PlanningRegion
% - RMRegion
% - SITEARCHITECTURE
% - Master
% - SPATIAL_REGION
% - MCC
% - MNC	

fid_ = fopen (fullfile(folder,site_atoll_filename),'r');
site_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
fclose(fid_);

%% Read LTE_BTS_trial_cluster_20100809
%The fields are:
% - Site
% - Transmitter
% - Active
% - DX (m)
% - DY (m)
% - Polarisation
% - Antenna	
% - Height (m)
% - Azimuth (°)
% - Mechanical Downtilt (°)
% - Main Calculation Radius (m)
% - Main Propagation Model
% - Transmission Loss (dB)
% - Reception Loss (dB)
% - Noise Figure (dB)
% - Hexagon groups
% - Hexagon radius (m)
% - Comments
% - TMA Equipment
% - Feeder Equipment
% - BTS Equipment
% - Transmission Feeder Length (m)
% - Reception Feeder Length (m)
% - Receiver antenna diversity gain (dB)
% - Miscellaneous Transmission Losses (dB)
% - Miscellaneous Reception Losses (dB)
% - Extended Calculation Radius (m)
% - Extended Propagation Model
% - Main Resolution (m)
% - Extended Resolution (m)
% - Additional Electrical Downtilt (°)
% - Number of Transmission Antenna Ports
% - Number of Reception Antenna Ports
% - Transmitter Type
% - COVERAGELAYER
% - BORDERCOORDINATION
% - Diversity
% - CoverageReported
% - PatternType
% - Optical	PARAMITRISATION_CLASS	

fid_ = fopen (fullfile(folder,transmitter_atoll_filename),'r');
transmitter_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
fclose(fid_);

%% Generate eNodeBs

% Now read the site and transmitter description and put it in the corresponding objects
% TODO : Change to print_log
% print_log(1,'Creating eNodeBs\n');
display('Creating eNodeBs');
eNodeBs = network_elements.eNodeB;

for site_idx = 1:size(site_atoll{1},1)
    
    id        = site_idx;
    name      = site_atoll{1}{site_idx};
    pos       = [ str2double(strrep(site_atoll{2}{site_idx},',','.')) str2double(strrep(site_atoll{3}{site_idx},',','.')) ]; % pos in meters (x,y)
    altitude  = str2double(site_atoll{4}{site_idx}(2:end-1));
    site_name = site_atoll{7}{site_idx};
    
    % Search for this eNodeB's sectors
    sectors_idx = search_string_in_textscan_output(transmitter_atoll{1},name);
    num_sectors = length(sectors_idx);
    
    % Fill in eNode object data
    eNodeBs(site_idx)           = network_elements.eNodeB;
    eNodeBs(site_idx).id        = id;
    eNodeBs(site_idx).name      = name;
    eNodeBs(site_idx).pos       = pos;
    eNodeBs(site_idx).altitude  = altitude;
    eNodeBs(site_idx).site_name = site_name;
    % neighbor info will be filled in a posteriori
    
    for sector_idx = 1:num_sectors
        if sector_idx==1
            eNodeBs(site_idx).sectors             = network_elements.eNodeB_sector;
        else
            eNodeBs(site_idx).sectors(sector_idx) = network_elements.eNodeB_sector;
        end
        
        transmitter_atoll_idx = sectors_idx(sector_idx);
        transmitter           = transmitter_atoll{2}{transmitter_atoll_idx};
        cell_idx              = find(strcmp(transmitter,cell_atoll{1}),1,'first');
        if isempty(cell_idx)
            error('The cell and transmitter Atoll files do not match. %s not found',transmitter);
        end
        switch capesso_params.planning_tool
            case 'atoll'
                % Do nothing
            case 'capesso'
                % Change the transmitter to read always the first sector
                % (we are using an isotrop antenna, so very probably we
                % will only have the pathloss map for one sector: no need
                % for the others).
                transmitter = regexprep(transmitter, '/U.*', '/U1');
            otherwise
                error('Only "atoll" and "capesso" supported');
        end

        antenna_name          = transmitter_atoll{7}{transmitter_atoll_idx};
        height                = str2double(strrep(transmitter_atoll{8}{transmitter_atoll_idx},',','.'));
        azimuth               = str2double(strrep(transmitter_atoll{9}{transmitter_atoll_idx},',','.'));
        %azimuth = 0;
        mechanical_downtilt   = str2double(strrep(transmitter_atoll{10}{transmitter_atoll_idx},',','.'));
        % Additional information contained in antenna_name
          antenna_name_split = textscan(antenna_name, '%s %s %f %f %f %s', 'delimiter', '|');
          antenna_type              = antenna_name_split{1}{1};
          % antenna_polarization    = antenna_name_split{2};
          % antenna_3dB_beamwidth_h = antenna_name_split{3};
          % antenna_3dB_beamwidth_v = antenna_name_split{4};
          electrical_downtilt       = antenna_name_split{5};
          % antenna_frequency_band  = antenna_name_split{6};
        % Kathrein antenna used for this sector
        antenna               = antennas.kathreinTSAntenna(kathrein_antenna_folder, antenna_type, frequency);
        power_dBm             = str2double(cell_atoll{3}{cell_idx});
        
        % Fill in eNodeB Sector Object data
        eNodeBs(site_idx).sectors(sector_idx).parent_eNodeB       = eNodeBs(site_idx);
        eNodeBs(site_idx).sectors(sector_idx).id                  = sector_idx;
        eNodeBs(site_idx).sectors(sector_idx).antenna             = antenna;
        eNodeBs(site_idx).sectors(sector_idx).max_power           = 10^(power_dBm/10) / 1000;
        eNodeBs(site_idx).sectors(sector_idx).transmitter         = transmitter;

        eNodeBs(site_idx).sectors(sector_idx).antenna_name        = antenna_name;
        eNodeBs(site_idx).sectors(sector_idx).antenna_type        = antenna_type;
        eNodeBs(site_idx).sectors(sector_idx).tx_height           = height;
        eNodeBs(site_idx).sectors(sector_idx).azimuth             = azimuth;
        eNodeBs(site_idx).sectors(sector_idx).mechanical_downtilt = mechanical_downtilt;        
        eNodeBs(site_idx).sectors(sector_idx).electrical_downtilt = electrical_downtilt;
    end
end

%% Some debug plotting
if capesso_params.enable_debug_plotting
    UEs = [];
    map_res = 20;
    current_TTI = 1;
    global LTE_config;
    LTE_config.plots.user_positions = 1;
    LTE_plot_show_network(eNodeBs, UEs, map_res, current_TTI);
end

function indexes = search_string_in_textscan_output(a_column,a_string)

indexes = [];

for i_=1:size(a_column,1)
    if strcmp(a_column{i_},a_string)
        indexes = [indexes i_];
    end
end
