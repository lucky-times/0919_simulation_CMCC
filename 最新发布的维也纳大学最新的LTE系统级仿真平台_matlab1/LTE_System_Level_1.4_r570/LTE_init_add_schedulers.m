function LTE_init_add_schedulers(eNodeBs,UEs,CQI_mapper,BLER_curves)
% Adds the needed scheduler type and resource block grid to each eNodeb's
% sector
% (c) Josep Colom Ikuno, INTHFT, 2008
% input:   eNodeBs  ... array of eNodeBs
%          UEs      ... array of UEs

global LTE_config;

switch LTE_config.scheduler
    case 'round robin'
        % Correct
    case 'best cqi'
        % Correct
    case 'proportional fair'
        % Correct
    case 'max min'
        % Correct
    case 'max TP'
        % Correct
    case 'resource fair'        
        % Correct
    case 'prop fair Sun'
        % Correct
    case 'constrained'
        % Correct
    case 'alpha fair'
        % Correct
    otherwise
        error([LTE_config.scheduler ' scheduler not supported']);
end

print_log(1,['Creating ' LTE_config.scheduler ' schedulers and resource block grids\n']);

% No reason to use a different SINR averager instance for each scheduler, we can reuse the same one
switch LTE_config.SINR_averaging.algorithm
    case 'MIESM'
        the_SINR_averager = utils.miesmAveragerFast(LTE_config.SINR_averaging.BICM_capacity_tables,LTE_config.SINR_averaging.betas);
    case 'EESM'
        the_SINR_averager = utils.eesmAverager(LTE_config.SINR_averaging.betas,LTE_config.SINR_averaging.MCSs);
    otherwise
        error('SINR averaging algorithm not supported');
end

% Add RB grid representation and scheduler to each sector.
% Set also homogeneous power load
for b_ = 1:length(eNodeBs)
    for s_=1:length(eNodeBs(b_).sectors)
        
        % Set whether the eNodeBs will always transmit, even if no UEs are attached.
        eNodeBs(b_).sectors(s_).always_on = LTE_config.always_on;
        max_data_power  = eNodeBs(b_).sectors(s_).max_power;
        signaling_power = eNodeBs(b_).sectors(s_).signaling_power;
        LTE_config.scheduler_params.max_power = max_data_power; % For backwards compatibility
        
        % Continue with Scheduler initialization
        switch LTE_config.scheduler
            case 'round robin'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.roundRobinScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
            case 'best cqi'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.bestCqiScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
            case 'max min'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.maxMinScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
            case 'max TP'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.maxThroughputScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
            case 'resource fair'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.resourceFairScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));    
            case 'prop fair Sun'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.propFairSunScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_)); 
            case 'constrained'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.ConstrainedScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
            case 'alpha fair'
                eNodeBs(b_).sectors(s_).scheduler = schedulers.alphaFairScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_)); 
            case 'proportional fair'
                if isfield(LTE_config,'scheduler_params') && isfield(LTE_config.scheduler_params,'alpha') && isfield(LTE_config.scheduler_params,'beta')
                    eNodeBs(b_).sectors(s_).scheduler = schedulers.proportionalFairScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_),LTE_config.scheduler_params.alpha,LTE_config.scheduler_params.beta);
                else
                    eNodeBs(b_).sectors(s_).scheduler = schedulers.proportionalFairScheduler(LTE_config.scheduler_params,eNodeBs(b_).sectors(s_));
                end
            otherwise
                error('Scheduler %s not defined',LTE_config.scheduler);
        end
        
        % Set scheduler SINR averaging algorithm
        eNodeBs(b_).sectors(s_).scheduler.SINR_averager = the_SINR_averager;

        % Other data required to perform SINR averaging at the transmitter side
        eNodeBs(b_).sectors(s_).scheduler.CQI_mapper = CQI_mapper;
        eNodeBs(b_).sectors(s_).scheduler.BLER_curves = BLER_curves;
        
        % Add genie information
        eNodeBs(b_).sectors(s_).scheduler.genie.UEs     = UEs;
        eNodeBs(b_).sectors(s_).scheduler.genie.eNodeBs = eNodeBs;
        
        % Add TTI delay information
        eNodeBs(b_).sectors(s_).scheduler.feedback_delay_TTIs = LTE_config.feedback_channel_delay;
        
        % TX modes:
        % 1: Single Antenna
        % 2: Transmit Diversity
        % 3: Open Loop Spatial Multiplexing
        % 4: Closed Loop SM
        
        % RB grid initialization
        eNodeBs(b_).sectors(s_).RB_grid = network_elements.resourceBlockGrid(LTE_config.N_RB,LTE_config.sym_per_RB_nosync,LTE_config.sym_per_RB_sync);
        eNodeBs(b_).sectors(s_).RB_grid.set_homogeneous_power_allocation(eNodeBs(b_).sectors(s_).max_power,eNodeBs(b_).sectors(s_).signaling_power);
        eNodeBs(b_).sectors(s_).RB_grid.tx_mode    = LTE_config.tx_mode;
    end
end

% Add each user to its corresponding scheduler (initialisation)
for u_=1:length(UEs)
    id = UEs(u_).id;
    b_ = UEs(u_).attached_eNodeB.id;
    s_ = UEs(u_).attached_sector;
    eNodeBs(b_).sectors(s_).scheduler.add_UE(id);
end
