function pregenerated_fast_fading = LTE_init_get_microscale_fading_SL_trace
% Generate the fading parameters that model the fast (microscale) fading at
% system level.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at.
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

%clc;
%close all;

%% Config

% Possible Tx modes are:
%   1: Single Antenna
%   2: Transmit Diversity
%   3: Open Loop Spatial Multiplexing
%   4: Closed Loop SM
% Number of antenna ports can either be 2 or 4
% Codebook index specifies the codebook index as of TS.36.211 (for closed loop SM)
% nLayers specifies how many layers (symbols) are going to be transmitted.
% Either 1, 2, 3 or 4

standalone = false;
show_trace_generation_debug = false;
global LTE_config;

if standalone
    LTE_config.debug_level            = 1;
    % Initial config (Winner II+)
    config.system_bandwidth           = 1.4e6;
    config.channel_type               = 'winner+';
    config.nTX                        = 2;
    config.nRX                        = 2;
    config.trace_length_s             = 1;
    config.UE_speed                   = 400/3.6; % converted to m/s. High speed ~ uncorrelated
    config.parallel_toolbox_installed = true; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = 0;
    config.correlated_fading          = true;
    config.f                          = 2.1400e+009;
    config.TTI_length                 = 1.0000e-003;
    config.tx_mode                    = 4; % OLSM
    
    winner_config                     = load('winner_trace_config');
    config.trace_params               = winner_config.config_trace_params;
    config.trace_params.speed         = config.UE_speed;
    clear winner_config
else
    % Initial config (simulator-linked)
    config.system_bandwidth           = LTE_config.bandwidth;
    config.channel_type               = LTE_config.channel_model.type;
    config.nTX                        = LTE_config.nTX;
    config.nRX                        = LTE_config.nRX;
    config.trace_length_s             = LTE_config.channel_model.trace_length;
    config.UE_speed                   = LTE_config.UE_speed; % converted to m/s
    config.parallel_toolbox_installed = LTE_config.parallel_toolbox_installed; % change it to false for testing purposes, as if not you will not be able to debug properly
    config.feedback_channel_delay     = LTE_config.feedback_channel_delay;
    config.correlated_fading          = LTE_config.channel_model.correlated_fading;
    config.f                          = LTE_config.frequency;
    config.trace_params               = LTE_config.trace_params;
    config.TTI_length                 = LTE_config.TTI_length;
    config.tx_mode                    = LTE_config.tx_mode;
end

% sigma_n2 = 10^((LTE_config.UE.receiver_noise_figure + LTE_config.UE.thermal_noise_density)/10)/1000;    % Receiver noise variance in Watt

% We now have all of the possible precoding combinations stored
precoding_configs  = get_all_precoding_combinations;
precoding_matrices = cell(1,length(precoding_configs));
for i_ = 1:length(precoding_configs)
    precoding_matrices{i_} = get_precoding_matrix(precoding_configs(i_));
end

% Channel trace for the target and interfering channels
switch config.channel_type
    case 'winner+'
        channel_factory_H0 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params);
        channel_factory_H1 = channel_gain_wrappers.winnerChannelFactory(config.system_bandwidth,config.trace_params);
    otherwise
        channel_factory_H0 = channel_gain_wrappers.pdpChannelFactory(config.system_bandwidth,config.trace_params);
        channel_factory_H1 = channel_gain_wrappers.pdpChannelFactory(config.system_bandwidth,config.trace_params);
end
print_log(1,sprintf('Generating %dx%d channel trace of length %3.2fs\\n',config.nTX,config.nRX,ceil(config.trace_length_s)));
H_trace0 = channel_factory_H0.generate_FF_trace(config.trace_length_s/config.TTI_length);

% Interfering channel trace
print_log(1,sprintf('Generating %dx%d interfering channel trace of length %3.2fs\\n',config.nTX,config.nRX,ceil(config.trace_length_s)));
H_trace1 = channel_factory_H1.generate_FF_trace(config.trace_length_s/config.TTI_length);

% Commented from LL vs SL
%H_trace0 = LTE_init_generate_FF_tracev2_LLvsSL(config.system_bandwidth,config.channel_type,config.nTX,config.nRX,config.trace_length_s,config.UE_speed,'UEchannel');
%H_trace1 = LTE_init_generate_FF_tracev2_LLvsSL(config.system_bandwidth,config.channel_type,config.nTX,config.nRX,config.trace_length_s,config.UE_speed,'interferers');

%% Channel normalization

% Note: each MIMO channel is normalized to a mean power of one
H_trace_normalized        = H_trace0.H_RB_samples;
H_trace_interf_normalized = H_trace1.H_RB_samples;

% Free up memory
clear H_trace0;
clear H_trace1;

% The channel trace is stored in H_trace.H_RB_samples, containing a
% (nRX,nTX,subframe_num,sample_num) matrix. We will calculate the SINR
% trace for each of these samples (every 6 subcarriers).

if config.parallel_toolbox_installed && ~matlabpool('size')
    matlabpool open;
end
switch config.tx_mode
    case 1
        % SISO trace
        SISO_trace = phy_modeling.txModeTrace;
        LTE_common_SISO_trace(config,H_trace_normalized,H_trace_interf_normalized,SISO_trace,show_trace_generation_debug);
    case 2
        % TxD (up to 2x2)
        clear precoding_matrix
        codebook_idx = 1; % 2 AT ports, 2 layers (in TxD nLayers = nAT ports)
        precoding_matrix = precoding_matrices{codebook_idx};
        
        TxD_2x2_trace = phy_modeling.txModeTrace;
        LTE_common_TxD_2x2_trace(config,H_trace_normalized,H_trace_interf_normalized,precoding_matrix,TxD_2x2_trace,show_trace_generation_debug);
    case 3
        % OLSM
        OLSM_trace = cell(min(config.nTX,config.nRX),1);
        switch config.nTX
            case 2
                sliced_precoding_matrices = precoding_matrices([3 4]);         % Precoding matrices for OLSM, 2 antenna ports and 1/2 layers
            otherwise
                error('2 TX antennas supported');
        end
        parfor layer_i = 1:min(config.nTX,config.nRX)
            OLSM_trace{layer_i} = LTE_common_OLSM_2x2_trace(config,H_trace_normalized,H_trace_interf_normalized,sliced_precoding_matrices{layer_i},show_trace_generation_debug,layer_i);
        end
    case 4
        % CLSM
        CLSM_trace = cell(min(config.nTX,config.nRX),1);
        switch config.nTX
            case 2
                sliced_precoding_matrices = precoding_matrices([8 9]);         % Precoding matrices for CLSM, 2 antenna ports and 1/2 layers
            case 4
                sliced_precoding_matrices = precoding_matrices([10 11 12 13]); % Precoding matrices for CLSM, 4 antenna ports and 1/2/3/4 layers
            otherwise
                error('2 or 4 TX antennas supported');
        end

        % Calculate the average H from each RB (needed for the PMI calculation)
        H_size   = size(H_trace_normalized); % channel matrix sizes
        N_sc     = 2;                        % number of subcarriers per resource block (just 2 samples are picked per RB)
        N_rb     = H_size(4)/N_sc;           % number of resource blocks
        H_t = zeros([H_size(1) H_size(2) H_size(3) N_rb]); % average channel of the RB
        
        for i_=1:H_size(3)
            H_TTI = reshape(H_trace_normalized(:,:,i_,:),[H_size(1) H_size(2) H_size(4)]);
            H_t(:,:,i_,:) = reshape(mean(reshape(H_TTI,[H_size(1)*H_size(2) N_sc N_rb]),2),[H_size(1) H_size(2) N_rb]);
        end
        
        %% Call trace generation function
        parfor layer_i = 1:min(config.nTX,config.nRX)
            CLSM_trace{layer_i} = LTE_common_CLSM_trace(config,H_trace_normalized,H_t,H_trace_interf_normalized,sliced_precoding_matrices{layer_i},show_trace_generation_debug);
        end
    otherwise
        error('Tx mode %d not supported',config.tx_mode);
end
if config.parallel_toolbox_installed
    matlabpool close;
end

%% Create output fast fading trace
pregenerated_fast_fading                      = phy_modeling.PregeneratedFastFading;
pregenerated_fast_fading.trace_length_s       = config.trace_length_s;
pregenerated_fast_fading.trace_length_samples = config.trace_length_s / 1e-3;
pregenerated_fast_fading.system_bandwidth     = config.system_bandwidth;
pregenerated_fast_fading.channel_type         = config.channel_type;
pregenerated_fast_fading.nTX                  = config.nTX;
pregenerated_fast_fading.nRX                  = config.nRX;
pregenerated_fast_fading.UE_speed             = config.UE_speed;

pregenerated_fast_fading.t_step               = 1e-3;
pregenerated_fast_fading.f_step               = 15e3*6;

switch config.tx_mode
    case 1
        % SISO trace (mode 1)
        pregenerated_fast_fading.traces{1} = SISO_trace;
    case 2
        % TxD trace (mode 2)
        pregenerated_fast_fading.traces{2} = TxD_2x2_trace;
    case 3
        % OLSM trace (mode 3)
        pregenerated_fast_fading.traces{3} = OLSM_trace;
    case 4
        % CLSM trace (mode 4)
        pregenerated_fast_fading.traces{4} = CLSM_trace;
end


function precoding_matrices = get_precoding_matrix(precoding_config)
% This function returns the precoding matrix for a specific transmission mode
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at, modified by Josep Colom
% Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

tx_mode = precoding_config.tx_mode;
nAtPort = precoding_config.nAtPort;
nLayers = precoding_config.nLayers;

%% Check correctness of the combination of nAtPort and nLayers
if tx_mode == 1
    % SISO case
    if (nAtPort~=1 && nLayers~=1)
        error('SISO mode only accepts 1 layers and 1 codeword');
    end
else
    if (nAtPort~=2 && nAtPort~=4)
        error('MIMO modes require of 2 o 4 antenna ports');
    end
    if nLayers > nAtPort
        error('number of layers must be equal or lower than number of antenna ports');
    end
end

if tx_mode == 2
    % TxD mode
    if (nAtPort ~= nLayers)
        error('TxD requires nLayers=nAtPorts');
    end
elseif tx_mode == 3
    % OLSM (named Large CDD in the standard)
elseif tx_mode == 4
    % CLSM
else
    error('tx_mode %d not defined in this function',tx_mode);
end

%% Codebook setting
% CLSM
if tx_mode == 4
    if ~isfield(precoding_config,'codebook_idxs')
        error('For tx_mode 4 a codebook index set must be specified');
    end
    codebook_indexs = precoding_config.codebook_idxs;
    % OLSM
elseif tx_mode == 3
    % OLSM uses codebooks 12-15 in a cyclic way. Thus we set codebook_index to [12 13 14 15]
    if nAtPort==4
        codebook_index = [12 13 14 15]-1;
    elseif nAtPort==2
        codebook_index = 0;
    end
end

LTE_params = LTE_params_function;

% Transmit diversity
if (tx_mode == 2)
    % We call the precoding matrix of TxD Z
    % Matrix corresponding to 36.211 section 6.3.4.3  
    precoding_matrices.Z = LTE_params.Z{nAtPort/2};
    precoding_matrices.name = 'TxD';
    precoding_matrices.tx_mode = tx_mode;
    precoding_matrices.nAtPort = nAtPort;
    precoding_matrices.nLayers = nLayers;

% Open loop spatial multiplexing, section 6.3.4.2.2 (Large CDD)
elseif (tx_mode == 3)   
    precoding_matrices.U = LTE_params.U_l{nLayers};
    precoding_matrices.D = LTE_params.D_l{nLayers};
    if (nAtPort == 2)
        W = LTE_params.W{nLayers}(:,:,codebook_index+1);
    else
        W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
        % nLayers long Cyclic precoding matrix
        W = zeros(4,nLayers,4);
        for ii = 13:16
            W(:,:,ii-12) = W_temp(:,LTE_params.mapping{nLayers}(ii,:),ii-12);
        end
    end
    precoding_matrices.W = W;
    precoding_matrices.name = 'OLSM';
    precoding_matrices.tx_mode = tx_mode;
    precoding_matrices.nAtPort = nAtPort;
    precoding_matrices.nLayers = nLayers;
    precoding_matrices.codebook_index = codebook_index;
    
% Closed loop spatial multiplexing, section 6.3.4.2.1
else
    W = zeros(nAtPort,nLayers,length(codebook_indexs));
    if (nAtPort == 2)
        if (min(codebook_indexs)<0 || max(codebook_indexs)>3) && nLayers ==1
            error('Only codebooks 0-3 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
        elseif (min(codebook_indexs)<0 || max(codebook_indexs)>2) && nLayers ==2
            error('Only codebooks 0-2 are defined for %d layers (see TS.36.211, Table 6.3.4.2.3-1)',nLayers);
        end
        for cb_ = 1:length(codebook_indexs)
            codebook_index = codebook_indexs(cb_);
            W(:,:,cb_) = LTE_params.W{nLayers}(:,:,codebook_index+1);
        end
    else
        if min(codebook_indexs)<0 || max(codebook_indexs)>15
            error('Only codebooks 0-15 are defined (see TS.36.211, Table 6.3.4.2.3-2)');
        end
        for cb_ = 1:length(codebook_indexs)
            codebook_index = codebook_indexs(cb_);
            W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
            W(:,:,cb_) = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
        end
    end
    precoding_matrices.W = W;
    precoding_matrices.name = 'CLSM';
    precoding_matrices.tx_mode = tx_mode;
    precoding_matrices.nAtPort = nAtPort;
    precoding_matrices.nLayers = nLayers;
    precoding_matrices.codebook_index = codebook_indexs;
end

function LTE_params = LTE_params_function
% Re-create needed load_parameters data from Link level for the generation of the precoding matrices.
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at


%% Create the Codebook for Precoding

% Transmit diversity
LTE_params.Z{1} =  [1, 0, 1i,  0;
         0,-1,  0, 1i;
         0, 1,  0, 1i;
         1, 0,-1i,  0];
LTE_params.Z{2} =  [1, 0, 0, 0, 1i,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0,-1, 0, 0,  0, 1i,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 1, 0, 0,  0, 1i,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         1, 0, 0, 0,-1i,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 1, 0,  0,  0, 1i, 0;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0,-1,  0,  0,  0,1i;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 0, 1,  0,  0,  0,1i;
         0, 0, 0, 0,  0,  0,  0, 0;
         0, 0, 1, 0,  0,  0,-1i, 0];
     
% Spatial multiplexing
U_temp = [  1,-1,-1,-1;     % Matrix corresponding to vectors u0 ... u15 in Table 6.3.4.2.3-2
            1,-1i,1,1i;
            1,1,-1,1;
            1,1i,1,-1i;
            1,(-1-1i)/sqrt(2), -1i,(1-1i)/sqrt(2);
            1,(1-1i)/sqrt(2), 1i,(-1-1i)/sqrt(2);
            1,(1+1i)/sqrt(2), -1i,(-1+1i)/sqrt(2);
            1,(-1+1i)/sqrt(2), 1i,(1+1i)/sqrt(2);
            1,-1,1,1;
            1,-1i,-1,-1i;
            1,1,1,-1;
            1,1i,-1,1i;
            1,-1,-1,1;
            1,-1,1,-1;
            1,1,-1,-1;
            1,1,1,1;].';
 Wn = zeros(4,4,16); 
 for ii = 1:16 
     LTE_params.Wn(:,:,ii)=diag(ones(1,4))-2*U_temp(:,ii)*U_temp(:,ii)'/(U_temp(:,ii)'*U_temp(:,ii));
 end
 
 % W Matrix according to Table 6.3.4.2.3-1
%  LTE_params.W{1} = cat(3,[1;0],[0;1],[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
%         [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
 LTE_params.W{1} = cat(3,[1/sqrt(2);1/sqrt(2)],[1/sqrt(2);-1/sqrt(2)],...
        [1/sqrt(2);1i/sqrt(2)],[1/sqrt(2);-1i/sqrt(2)]);
 LTE_params.W{2} = cat(3,1/sqrt(2)*[1,0;0,1],1/(2)*[1,1;1,-1],1/(2)*[1,1;1i,-1i]);
 
 % Large delay CDD  
 LTE_params.U_l{1} = 1;
 LTE_params.U_l{2} = 1/sqrt(2)*[1,1;1,exp(-1i*pi)]; 
 LTE_params.U_l{3} = 1/sqrt(3)*[1,1,1;1,exp(-1i*2*pi/3),exp(-1i*4*pi/3);1,exp(-1i*4*pi/3),exp(-1i*8*pi/3)];
 LTE_params.U_l{4} = 1/2*[1,1,1,1;1,exp(-1i*2*pi/4),exp(-1i*4*pi/4),exp(-1i*6*pi/4);...
                            1,exp(-1i*4*pi/4),exp(-1i*8*pi/4),exp(-1i*12*pi/4);...
                            1,exp(-1i*6*pi/4),exp(-1i*12*pi/4),exp(-1i*18*pi/4)];
 LTE_params.D_l{1} = 1;
 LTE_params.D_l{2} = [1,0;0,exp(-1i*pi)];
 LTE_params.D_l{3} = [1,0,0;0,exp(-1i*2*pi/3),0;0,0,exp(-1i*4*pi/3)];
 LTE_params.D_l{4} = [1,0,0,0;0,exp(-1i*2*pi/4),0,0;0,0,exp(-1i*4*pi/4),0;0,0,0,exp(-1i*6*pi/4)];
 
 % Note that as of v.8.3.0, small delay CDD is removed from the standard
 % (28/05/08	RAN_40	RP-080432	0043	-	Removal of small-delay CDD
 
 % Precoding matrix W columns to take for each layer mapping
 LTE_params.mapping{1} = ones(16,1);
 LTE_params.mapping{2}=[1 4;1 2;1 2;1 2;1 4;1 4;1 3;1 3;1 2;1 4;1 3;1 3;1 2;1 3;1 3;1 2];
 LTE_params.mapping{3}=[1 2 4;1 2 3;1 2 3;1 2 3;1 2 4;1 2 4;1 3 4;1 3 4;1 2 4;1 3 4;1 2 3;1 3 4;1 2 3;1 2 3;1 2 3;1 2 3];
 LTE_params.mapping{4}=[1 2 3 4;1 2 3 4;3 2 1 4;3 2 1 4;1 2 3 4;1 2 3 4;1 3 2 4;
     1 3 2 4;1 2 3 4;1 2 3 4;1 3 2 4;1 3 2 4;1 2 3 4;1 3 2 4;3 2 1 4;1 2 3 4];
 
function precoding_config = get_all_precoding_combinations
% This small helper function returns all possible precoding options for LTE.
% (c) Josep Colom Ikuno, INTHFT, 2008
% www.nt.tuwien.ac.at

% TxD
precoding_config.tx_mode = 2;
precoding_config.nAtPort = 2;
precoding_config.nLayers = 2;

precoding_config(2).tx_mode = 2;
precoding_config(2).nAtPort = 4;
precoding_config(2).nLayers = 4;

% Large delay CDD (OLSM)

% No mode with one layer/rank 1 (TxD is used in that case)
precoding_config(3).tx_mode = 3;
precoding_config(3).nAtPort = 2;
precoding_config(3).nLayers = 1;

precoding_config(4).tx_mode = 3;
precoding_config(4).nAtPort = 2;
precoding_config(4).nLayers = 2;

precoding_config(5).tx_mode = 3;
precoding_config(5).nAtPort = 4;
precoding_config(5).nLayers = 2;

precoding_config(6).tx_mode = 3;
precoding_config(6).nAtPort = 4;
precoding_config(6).nLayers = 3;

precoding_config(7).tx_mode = 3;
precoding_config(7).nAtPort = 4;
precoding_config(7).nLayers = 4;

% CLSM
precoding_config(8).tx_mode = 4;
precoding_config(8).nAtPort = 2;
precoding_config(8).nLayers = 1;
precoding_config(8).codebook_idxs = 0:3;

precoding_config(9).tx_mode = 4;
precoding_config(9).nAtPort = 2;
precoding_config(9).nLayers = 2;
precoding_config(9).codebook_idxs = 1:2;


precoding_config(10).tx_mode = 4;
precoding_config(10).nAtPort = 4;
precoding_config(10).nLayers = 1;
precoding_config(10).codebook_idxs = 0:15;

precoding_config(11).tx_mode = 4;
precoding_config(11).nAtPort = 4;
precoding_config(11).nLayers = 2;
precoding_config(11).codebook_idxs = 0:15;

precoding_config(12).tx_mode = 4;
precoding_config(12).nAtPort = 4;
precoding_config(12).nLayers = 3;
precoding_config(12).codebook_idxs = 0:15;

precoding_config(13).tx_mode = 4;
precoding_config(13).nAtPort = 4;
precoding_config(13).nLayers = 4;
precoding_config(13).codebook_idxs = 0:15;
