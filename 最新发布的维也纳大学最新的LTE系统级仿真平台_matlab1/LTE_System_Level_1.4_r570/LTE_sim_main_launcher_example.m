close all force;
clc;
%clear all
%clear global;
%clear classes;

%% Load parameters. Now done outside
LTE_load_params_hex_grid_tilted;

%% If you want to modify something taking as a base the configuration file, do it here: here an example is show that changes the inter-eNodeB distances based on the LTE_load_params_hex_grid_tilted config file.
%LTE_config.cache_network = false;
%LTE_config.network_cache = 'auto';
%LTE_config.inter_eNodeB_distance = 500;
%LTE_config.map_resolution = 20;
%LTE_config.simulation_time_tti = 25; % simulate 25 TTIs
%LTE_config.results_file           = 'auto';

% LTE_load_params_dependant;

print_log(1,'Loaded configuration file\n');
LTE_sim_main
