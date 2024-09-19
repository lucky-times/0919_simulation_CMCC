function [PMI] = LTE_feedback(H,nLayers,W_all,nAtPort,receiver) 
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the PMI

%% some initializations
N_sc = 2;                   % number of subcarriers per resource block (just 2 samples are picked per RB)
N_rb = size(H,4)/N_sc;      % number of resource blocks  
PMI = zeros(N_rb,1);        % number of PMI values    
switch nLayers              % precoder indices
    case 1
        if nAtPort == 2
           i_loop = 0:3;
           add = 1;
        else 
           i_loop = 0:15;
           add = 1;
        end
    case 2
        if nAtPort == 2
           i_loop = 1:2;
           add = 0;
        else 
           i_loop = 0:15;
           add = 1;
        end
    case {3,4}
        i_loop = 0:15;
        add = 1;
end

for RB_i = 1:N_rb
        freq_band = (RB_i-1)*N_sc+1:RB_i*N_sc;
        H_t = mean(H(:,:,1,freq_band),4);   % average channel of current RB
        I = zeros(length(i_loop),1);   
        SNR = zeros(length(i_loop),nLayers);
        for i = i_loop
            W = W_all(:,:,i+add);
            P = H_t*W;
            switch receiver     % right now only ZF receiver supported in SL
                case 'ZF'
                    F = pinv(P);   % ZF receiver                  
                case 'MMSE'
                    temp = P'*P;
                    F = (temp+0.01*eye(size(temp)))^-1*P';  % MMSE receiver
                otherwise
                    temp = P'*P;
                    F = (temp+0.01*eye(size(temp)))^-1*P';  % MMSE receiver
            end
            K = F*P;
            SNR(i+add,1:nLayers) = abs(diag(K)).^2./(sum(abs(K-diag(diag(K))).^2,2)+0.01.*sum(abs(F).^2,2));
            I(i+add) = sum(log2(1+SNR(i+add,:)));  % rate of one resource block for different precoders and ranks
        end
    [~,C1] = max(I,[],1);
    PMI(RB_i) = C1;
end

