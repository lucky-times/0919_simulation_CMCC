close all force;
clc;
% clear global;
% clear classes;
cd ..

for loop = 1:5
    clearvars -except calling_dir SL_dir LL_dir loop
    switch loop
        case 1
            scheduler = 'best cqi';
            output_filename = 'SL_BCQI';
        case 2
            scheduler = 'max min';
            output_filename = 'SL_MM';
        case 3
            scheduler = 'round robin';
            output_filename = 'SL_RR';
        case 4
            scheduler = 'prop fair Sun';
            output_filename = 'SL_PF';
        case 5
            scheduler = 'resource fair';
            output_filename = 'SL_RF';
    end
    LTE_load_params_scheduler_comp;
    print_log(1,'Loaded configuration file\n');
    LTE_sim_main
    save([output_filename '.mat'],'simulation_traces');
end