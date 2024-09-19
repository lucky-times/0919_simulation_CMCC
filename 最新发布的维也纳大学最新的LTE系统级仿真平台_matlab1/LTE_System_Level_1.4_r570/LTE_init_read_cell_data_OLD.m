function eNodeBs = LTE_init_read_atoll_cell_data(capesso_params)
% Reads out the specified folder and searchs for the following files:
%   - *.par files 
%   - transmitter_atoll.txt
%   - Site_atoll.txt
%   - cell_atoll.txt
%
% (c) Josep Colom Ikuno, Martin Taranetz, INTHFT, 2010

%% Parameter definitions
folder                        = capesso_params.pathloss_data_folder;
cell_atoll_filename           = capesso_params.cell_atoll_filename;
site_atoll_filename           = capesso_params.site_atoll_filename;
transmitter_atoll_filename    = capesso_params.transmitter_atoll_filename;
kathrein_antenna_folder       = capesso_params.kathrein_antenna_folder;

% Workaround for using frequency defined in LTE_load_params
if capesso_params.frequency == 2e+9
    frequency = 2140;
else
    error('No antenna for demanded frequency');
end

%% Read out cell_atoll

% The fields are
%   - Name
%   - Transmitter
%   - Carrier
%   - Scrambling Code Domain
%   - Primary Scrambling Code
%   - AS Threshold (dB)
%   - Max Power (dBm)
%   - Pilot Power (dBm)
%   - SCH/Pilot Offset(dB)
%   - Other CCH/Pilot Offset (dB)
%   - Total Power (dBm)
%   - UL Load Factor (%)
%   - Comments
%   - SC Reuse Distance (m)
%   - Max Number of Intra-carrier Neighbours
%   - Max number of inter-technology neighbours
%   - Max UL Load Factor (%)
%   - Max DL Load (% Pmax)
%   - HSPA
%   - Available HSDPA Power (dBm)
%   - Power Headroom (dB)
%   - Max Number of HS-PDSCH Codes
%   - Max Number of Inter-carrier Neighbours
%   - Min Number of HS-PDSCH Codes
%   - HSDPA Dynamic Power Allocation
%   - HS-SCCH Dynamic Power Allocation
%   - HS-SCCH/Pilot Offset (dB)
%   - Number of HS-SCCH Channels
%   - Max Number of HSDPA Users
%   - HSDPA Scheduler Algorithm
%   - MUG Table=f(No. Users)
%   - UL Peak Rate per User (kbps)
%   - DL Peak Rate per User (kbps)
%   - Active
%   - HSUPA
%   - DL HSUPA/Pilot Offset (dB)
%   - Max Number of HSUPA Users
%   - UL Load Factor due to HSUPA (%)
%   - Number of HSUPA Users
%   - Reuse Factor (UL)
%   - Number of HSDPA Users

% Since the numerical fields are separated by a COMMA instead of a point, I
% will have to do the 4,7 -> 4.7 conversion by hand
fid_ = fopen (fullfile(folder,cell_atoll_filename),'r');
% M = textscan(fid_,'%s %s %n %s %n %n %n %n %n %s %n %n %s %n %n %n %n %n %s %s %n %n %n %n %s %s %n %n %n %s %s %n %n %s %s %n %n %n %n %n %n','Delimiter','\t','Headerlines',1);
cell_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
fclose(fid_);

%% Read site_atoll

% The fields are:
%   - Name
%   - X
%   - Y
%   - Altitude (m)
%   - Comments
%   - Max No. of UL CEs
%   - Max No. of DL CEs
%   - Equipment
%   - PYLON_HEIGHT (m)
%   - SUPPORT_NATURE
%   - SITENAME
%   - PlanningRegion
%   - RMRegion
%   - SITEARCHITECTURE
%   - VENDOR
%   - SITESTATUS
%   - Master
%   - CE_R99
%   - CE_HSDPA
%   - CE_HSUPA
%   - SPATIAL_REGION
%   - MCC
%   - MNC

fid_ = fopen (fullfile(folder,site_atoll_filename),'r');
site_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
fclose(fid_);

%% Read transmitter_atoll

%The fields are:
%   - Site
%   - Transmitter
%   - Active
%   - DX (m)
%   - DY (m)
%   - Frequency Band
%   - Antenna
%   - Height (m)
%   - Azimuth (°)
%   - Mechanical Downtilt (°)
%   - Main Calculation Radius (m)
%   - Main Propagation Model
%   - Transmission Loss (dB)
%   - Reception Loss (dB)
%   - Noise Figure (dB)
%   - Hexagon Groups
%   - Hexagon Radius (m)
%   - Comments
%   - TMA Equipment
%   - Feeder Equipment
%   - BTS Equipment
%   - Transmission Feeder Length (m)
%   - Reception Feeder Length (m)
%   - Receiver Antenna Diversity Gain (dB)
%   - Miscellaneous Transmission Losses (dB)
%   - Miscellaneous Reception Losses (dB)
%   - Extended Calculation Radius (m)
%   - Extended Propagation Model
%   - Main Resolution (m)
%   - Extended Resolution (m)
%   - Cell Edge Coverage Probability (%)
%   - Additional Electrical Downtilt (°)
%   - Reception Diversity
%   - Transmission Diversity
%   - Inter-carrier Power Sharing
%   - Max Shared Power (dBm)
%   - COVERAGELAYER
%   - BORDERCOORDINATION
%   - Diversity
%   - CoverageReported
%   - PatternType
%   - Optical
%   - PARAMITRISATION_CLASS

fid_ = fopen (fullfile(folder,transmitter_atoll_filename),'r');
transmitter_atoll = textscan(fid_,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t','Headerlines',1);
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
    site_name = site_atoll{11}{site_idx};
    
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
        frequency_band        = transmitter_atoll{6}{transmitter_atoll_idx};
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
        
        % Fill in eNodeB Sector Object data
        eNodeBs(site_idx).sectors(sector_idx).parent_eNodeB       = eNodeBs(site_idx);
        eNodeBs(site_idx).sectors(sector_idx).id                  = sector_idx;
        eNodeBs(site_idx).sectors(sector_idx).antenna             = antenna;
        eNodeBs(site_idx).sectors(sector_idx).max_power           = capesso_params.eNodeB_sector_tx_power;
        eNodeBs(site_idx).sectors(sector_idx).transmitter         = transmitter;
        eNodeBs(site_idx).sectors(sector_idx).frequency_band      = frequency_band;
        eNodeBs(site_idx).sectors(sector_idx).antenna_name        = antenna_name;
        eNodeBs(site_idx).sectors(sector_idx).antenna_type        = antenna_type;
        eNodeBs(site_idx).sectors(sector_idx).height              = height;
        eNodeBs(site_idx).sectors(sector_idx).azimuth             = azimuth;
        eNodeBs(site_idx).sectors(sector_idx).mechanical_downtilt = mechanical_downtilt;        
        eNodeBs(site_idx).sectors(sector_idx).electrical_downtilt = electrical_downtilt;
    end
end

%% A posteriori calculation of neighboring eNodeBs
for b_ = 1:length(eNodeBs)
    eNodeBs(b_).neighbors = LTE_init_get_eNodeB_neighbors(eNodeBs(b_), eNodeBs , capesso_params.inter_bts_distance+1);
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
