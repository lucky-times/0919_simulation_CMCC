close all force;
clc;
% clear global;
% clear classes;
cd ..

for loop = 1:4
    clearvars -except calling_dir SL_dir LL_dir loop
    switch loop
        case 1 % SISO
            T_TIs = 10000;
            mode = 1;
            nTX = 1;
            nRX = 1;
            position_filename = 'SLvsLL_SISO_comparison.mat';
            channel_filename = 'SLvsLL_SISO_TU.mat';
            output_filename = 'SISO_TU_SL';
        case 2 % TxD
            T_TIs = 16500;
            mode = 2;
            nTX = 2;
            nRX = 2;
            position_filename = 'SLvsLL_TxD_comparison.mat';
            channel_filename = 'SLvsLL_TxD_TU.mat';
            output_filename = 'TxD_TU_SL';
        case 3 % OLSM
            T_TIs = 33000;
            mode = 3;
            nTX = 2;
            nRX = 2;
            position_filename = 'SLvsLL_OLSM_comparison.mat';
            channel_filename = 'SLvsLL_OLSM_TU.mat';
            output_filename = 'OLSM_TU_SL';
        case 4 % CLSM
            T_TIs = 20000;
            mode = 4;
            nTX = 4;
            nRX = 2;
            position_filename = 'SLvsLL_CLSM_comparison.mat';
            channel_filename = 'SLvsLL_CLSM_TU.mat';
            output_filename = 'CLSM_TU_SL';
    end
    LTE_load_params_LLvsSL;
    print_log(1,'Loaded configuration file\n');
    LTE_sim_main
    save([output_filename '.mat'],'simulation_traces');
end