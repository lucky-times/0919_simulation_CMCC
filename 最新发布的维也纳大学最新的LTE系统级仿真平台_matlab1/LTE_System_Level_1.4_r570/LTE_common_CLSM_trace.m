function [ trace_to_fill ] = LTE_common_CLSM_trace( config,H_trace_normalized,H_t,H_trace_interf_normalized,precoding_matrix,debug_output)
% Generate the fading trace for the 2x2 CLSM LTE mode using the PMI
% feedback proposed in "Mutual Information based calculation of the
% precoding matrix indicator for 3GPP UMTS/LTE", ITG WSA 2010
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at.
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

% Re-create config params from input
system_bandwidth = config.system_bandwidth;
channel_type     = config.channel_type;
nTX              = config.nTX;
nRX              = config.nRX;
trace_length_s   = config.trace_length_s;

fprintf('CLSM %dx%d: %d Layers\n',nTX,nRX,precoding_matrix.nLayers);

UE_speed         = config.UE_speed;
nLayers          = precoding_matrix.nLayers;

% Preallocate
TTI_length = trace_length_s*1000;
SC_samples = size(H_trace_normalized,4);

H0F       = zeros(nRX,nLayers,SC_samples,TTI_length);
H1F       = zeros(nRX,nLayers,SC_samples,TTI_length);
H0F_pinv  = zeros(nLayers,nRX,SC_samples,TTI_length);
A         = zeros(nLayers,nLayers,SC_samples,TTI_length);
C         = zeros(nLayers,nLayers,SC_samples,TTI_length);
PMI_trace = zeros(SC_samples/2,TTI_length);
PMI       = zeros(SC_samples/2,TTI_length+config.feedback_channel_delay);

for TTI_ = 1:TTI_length+config.feedback_channel_delay
    % This construction allows a parfor construct to be used here
    if TTI_ <= TTI_length
        TTI_H         = H_trace_normalized(:,:,TTI_,:);        
        TTI_H_inter   = H_trace_interf_normalized(:,:,TTI_,:);
        TTI_H_t       = reshape(H_t(:,:,TTI_,:),[size(H_t,1) size(H_t,2) size(H_t,4)]);
    else
        TTI_H         = H_trace_normalized(:,:,mod(TTI_,TTI_length),:);   
        TTI_H_inter   = H_trace_interf_normalized(:,:,mod(TTI_,TTI_length),:);
        TTI_H_t       = reshape(H_t(:,:,mod(TTI_,TTI_length),:),[size(H_t,1) size(H_t,2) size(H_t,4)]);
    end
    [PMI(:,TTI_)] = LTE_feedback(TTI_H_t,nLayers,precoding_matrix.W,nTX,'ZF');
    
    if TTI_ >  config.feedback_channel_delay   
        TTI__ = TTI_ - config.feedback_channel_delay;
        PMI_act = PMI(:,TTI__);

        [   H0F(:,:,:,TTI__)...
            H1F(:,:,:,TTI__)...
            H0F_pinv(:,:,:,TTI__)...
            A(:,:,:,TTI__)...
            C(:,:,:,TTI__)...
            PMI_trace(:,TTI__)] = calculate_TTI_params(TTI_H,TTI_H_inter,precoding_matrix.W,nRX,nLayers,SC_samples,PMI_act);
    end
end


%% Extract now the fading parameters
zeta  = zeros(nLayers,SC_samples,TTI_length);  % Scales the received signal (1 for perfect channel knowledge)
for layer_idx = 1:nLayers
    zeta(layer_idx,:,:) = squeeze(abs(A(layer_idx,layer_idx,:,:)).^2);
end
chi   = reshape(sum(abs(A).^2,2),size(zeta,1),size(zeta,2),size(zeta,3)) - zeta;             % Represents inter-layer interference (0 for perfect channel knowledge)
psi   = reshape(sum(abs(H0F_pinv).^2,2),size(H0F_pinv,1),size(H0F_pinv,3),size(H0F_pinv,4)); % Scales the noise
theta = reshape(sum(abs(C).^2,2),size(C,1),size(C,3),size(C,4));                             % Scales the interference


%% Some testing (specially norms)

if debug_output
    norms_A = zeros(1,TTI_length);
    norms_H0_normalized = zeros(1,TTI_length);
    norms_H1_normalized = zeros(1,TTI_length);
    for TTI_ = 1:TTI_length
        norms_H0(TTI_) = norm(H_trace_normalized(:,:,TTI_,1),'fro')^2;
        norms_H1(TTI_) = norm(H_trace_interf_normalized(:,:,TTI_,1),'fro')^2;
        norms_H0F(TTI_)= norm(H0F(:,:,1,TTI_),'fro')^2;
        norms_A(TTI_)  = norm(A(:,:,1,TTI_),'fro')^2;
        norms_B(TTI_)  = norm(H0F_pinv(:,:,1,TTI_),'fro')^2;
        norms_C(TTI_)  = norm(C(:,:,1,TTI_),'fro')^2;
    end
    fprintf('<||H0||^2> = %3.2f\n',mean(norms_H0));
    fprintf('<||H1||^2> = %3.2f\n',mean(norms_H1));
    fprintf('<||F||^2>  = %3.2f\n',norm(precoding_matrix.W(:,:,1),'fro')^2);  % NOTE: just any precoder
    fprintf('<||H0F||^2>  = %3.2f\n',mean(norms_H0F));
    fprintf('<||A||^2>  = %3.2f\n',mean(norms_A));
    fprintf('<||B||^2>  = %3.2f\n',mean(norms_B));
    fprintf('<||C||^2>  = %3.2f\n',mean(norms_C));
    
    fprintf('Averages:\n');
    fprintf(' psi:   %3.2f\n',mean(psi(:)));
    fprintf(' theta: %3.2f\n',mean(theta(:)));
end

%% Fill in the output trace object
trace_to_fill                  = phy_modeling.txModeTrace;
trace_to_fill.tx_mode          = 4;
trace_to_fill.trace_length_s   = trace_length_s;
trace_to_fill.system_bandwidth = system_bandwidth;
trace_to_fill.channel_type     = channel_type;
trace_to_fill.nTX              = nTX;
trace_to_fill.nRX              = nRX;
trace_to_fill.UE_speed         = UE_speed;

trace_to_fill.trace.zeta  = zeta;
trace_to_fill.trace.chi   = chi;
trace_to_fill.trace.psi   = psi;
trace_to_fill.trace.theta = theta;
trace_to_fill.trace.PMI = PMI_trace;

%% Some plotting

if debug_output
    for i_=1:nLayers
        figure;
        hold on;
        plot(norms_H0,'k','Displayname','||H_0||_{F}^2 (Channel norm)');
        plot(norms_H1,'k:','Displayname','||H_1||_{F}^2 (Interf channel norm)');
        plot(squeeze(zeta(i_,1,:)),'r','Displayname','\zeta (RX power)');
        plot(squeeze(chi(i_,1,:)),'g','Displayname','\chi (inter-layer interf)');
        plot(squeeze(psi(i_,1,:)),'b','Displayname','\psi (noise enhancement)');
        plot(squeeze(theta(i_,1,:)),'m','Displayname','\theta (inter-cell interference)');
        set(gca,'Yscale','log');
        title(sprintf('CLSM, Layer %d/%d, subcarrier 1',i_,nLayers));
        grid on;
        legend('show','Location','best');
    end
end

function [H0F H1F H0F_pinv A C PMI] = calculate_TTI_params(H_trace_normalized,H_trace_interf_normalized,W,nRX,nLayers,SC_samples,PMI)
% H_trace_normalized(:,:,1,N_SCs)
% H_trace_interf_normalized(:,:,1,N_SCs)

H0F      = zeros(nRX,nLayers,1,SC_samples);
H1F      = zeros(nRX,nLayers,1,SC_samples);
H0F_pinv = zeros(nLayers,nRX,1,SC_samples);
A        = zeros(nLayers,nLayers,1,SC_samples);
C        = zeros(nLayers,size(H1F,2),1,SC_samples);

for SC_sample = 1:SC_samples
    current_H0 = H_trace_normalized(:,:,1,SC_sample);
    current_H1 = H_trace_interf_normalized(:,:,1,SC_sample);
    W_idx_H0 = PMI(ceil(SC_sample/2)); % H0 precoding matrix
    W_idx_H1 = 1;                      % dummy interfering precoding matrix
    
    % Get the pinv(HW) matrix for each TTI, SC sample and precoding matrix
    H0F_      = current_H0 * W(:,:,W_idx_H0);
    H1F_      = current_H1 * W(:,:,W_idx_H1);
    H0F_pinv_ = pinv(H0F_);
    
    H0F(:,:,1,SC_sample)      = H0F_;
    H1F(:,:,1,SC_sample)      = H1F_;
    
    H0F_pinv(:,:,1,SC_sample) = pinv(H0F(:,:,1,SC_sample));
    A(:,:,1,SC_sample)        = H0F_pinv_ * H0F_;
    C(:,:,1,SC_sample)        = H0F_pinv_ * H1F_;
end