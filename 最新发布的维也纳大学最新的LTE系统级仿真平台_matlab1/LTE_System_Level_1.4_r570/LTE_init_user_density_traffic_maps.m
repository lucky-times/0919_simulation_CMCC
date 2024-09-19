function user_density_traffic_maps = LTE_init_user_density_traffic_maps(udtm_folder, udtm_file_name, udtm_environment, enable_plotting)
% Martin Taranetz, INTHFT 2010

% Function returns user density traffic maps and descriptive data
% 4 traffic maps (environments) with different penetration losses are available:
%                           Penetration Loss[dB]
% o)  DI    Deep Indoor     23
% o)  ID    Indoor          17
% o)  IC    Incar           7
% o)  OD    Outdoor         0

global LTE_config;


%% Read out UDTM from .bil file and additional information from .HDR file
% udtm ... user traffic density map
filename_ = strcat(udtm_file_name,'_',udtm_environment,'_');
data_file = strcat(filename_,'focus.bil');
hdr_file  = strcat(filename_,'focus.HDR');
udtm_file_hdr = fopen(fullfile(udtm_folder, hdr_file) ,'r');
udtm_file     = fopen(fullfile(udtm_folder, data_file),'r');

% Header file information
udtm_hdr = textscan(udtm_file_hdr, '%s %f');

udtm.description.NWxmap = udtm_hdr{2}(1);
udtm.description.NWymap = udtm_hdr{2}(2);
udtm.description.xdim   = udtm_hdr{2}(3);
udtm.description.ydim   = udtm_hdr{2}(4);
udtm.description.ncols  = udtm_hdr{2}(5);
udtm.description.nrows  = udtm_hdr{2}(6);
udtm.description.nbits  = udtm_hdr{2}(7);
udtm.description.nbands = udtm_hdr{2}(8);

udtm.description.SWxmap = udtm.description.NWxmap;                                                % Lower-leftmost corner (x) -> SW
udtm.description.SWymap = udtm.description.NWymap - udtm.description.ydim*(udtm.description.nrows-1/100); % Lower-leftmost corner (y) -> SW
   
udtm.description.NExmap = udtm.description.NWxmap + udtm.description.xdim*(udtm.description.ncols-1/100);
udtm.description.NEymap = udtm.description.NWymap;
    
udtm.description.SExmap = udtm.description.NExmap;
udtm.description.SEymap = udtm.description.SWymap;
    
udtm.description.roi_x = [udtm.description.SWxmap udtm.description.SExmap];
udtm.description.roi_y = [udtm.description.SWymap udtm.description.NWymap];

% Read the .bil file with the appropriate coding
if udtm.description.nbits == 32
   udtm_map_t = fread(udtm_file, [udtm.description.ncols, udtm.description.nrows], 'single','l'); % Directly specify to use little endian, just in case...
   udtm_map_t(udtm_map_t==-realmax('single')) = 0; % In order to fill the 'holes' in the traffic map
else
   error('Precision of User Traffic Denisty Map not supported');
end

udtm_map_t     = udtm_map_t*LTE_config.traffic_map_upscaling;
udtm.data      = flipud(udtm_map_t');
udtm.upscaling = LTE_config.traffic_map_upscaling;

% Plotting
if enable_plotting
   figure;
   imagesc(udtm.description.roi_x,udtm.description.roi_y, udtm.data);
   hold on;
   set(gca,'YDir','normal');
   xlabel('x pos [m]');
   ylabel('y pos [m]');
   colorbar;
   caxis([0 10]); % Using 'jet' colormap this fits perfectly to MKA legend
   title(['User density traffic map (m) [Users/km^2] ' udtm_environment]);
end

fclose(udtm_file);
fclose(udtm_file_hdr);

user_density_traffic_maps = udtm;
