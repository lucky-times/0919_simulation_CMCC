close all force;
clc;
clear
clear global;
clear classes;

LTE_load_params_example_capesso;
print_log(1,'Loaded configuration file\n');
LTE_sim_main