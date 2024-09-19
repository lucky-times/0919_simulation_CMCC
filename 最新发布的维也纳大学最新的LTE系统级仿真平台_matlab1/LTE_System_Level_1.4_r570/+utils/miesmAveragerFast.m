classdef miesmAveragerFast < utils.sinrAverager
    % Implements a MIESM averager. Needs a precalculated table containing
    % the BICM capacity for each of the constellations used, which are here
    % classified according to their modulation orders. Alternative
    % implementation to the original MIESM SINR averager that works faster
    % This MIESM averager DOES support calibration
    % (c) Josep Colom Ikuno, INTHFT, 2010

   properties
       length
       BICM_matrix_x
       BICM_matrix_y
       BICM_matrix_inv_x
       BICM_matrix_inv_y
       resolution = 0.01;
       resolution_inv
       multipliers
       SINR_range
       BICM_range
       nCQIs
       betas
       betas_dB
   end

   methods
       function obj = miesmAveragerFast(mat_file_to_load,betas,varargin)
           % Create the inverse I curve (the same as with the previous implementation)
           load(mat_file_to_load);
           
           % Upsample the BICM curves
           SNR_vector = BICM_capacity_tables(1).SNR(1):obj.resolution:BICM_capacity_tables(1).SNR(end);
           
           CQI_range = LTE_common_get_CQI_params('range');
           CQI_tables = LTE_common_get_CQI_params(CQI_range(1):CQI_range(2));
           n_CQIs = length(CQI_tables);
           mod_orders = [CQI_tables.modulation_order];
           for i_=1:length(BICM_capacity_tables)
               % Overwrite loaded values with upsampled ones
               BICM_capacity_tables(i_).I   = interp1(BICM_capacity_tables(i_).SNR,BICM_capacity_tables(i_).I,SNR_vector);
               BICM_capacity_tables(i_).SNR = SNR_vector;
               
               % Generate the inverse mapping
               BICM_capacity_tables(i_).I(1) = 0;
               BICM_capacity_tables(i_).I(end) = ceil(BICM_capacity_tables(i_).I(end));
               [b,m,n]=unique(BICM_capacity_tables(i_).I,'first');
               
               % To have only unique values
               BICM_capacity_tables(i_).I_inv = b;
               BICM_capacity_tables(i_).SNR_inv = BICM_capacity_tables(i_).SNR(m);
           end
           
           % Assume that the SNR is the same one for all data. It should be
           % automaticallythe case if the BICM capacity script of this
           % simulator was used.
           min_C = 0;
           max_C = max([BICM_capacity_tables.I]);
           BICM_matrix_x     = SNR_vector;
           BICM_matrix_y     = zeros(length(BICM_matrix_x),n_CQIs);
           BICM_matrix_inv_x = linspace(min_C,max_C,length(BICM_matrix_x));
           BICM_matrix_inv_y = zeros(length(BICM_matrix_inv_x),n_CQIs);

           obj.nCQIs = length(CQI_range(1):CQI_range(2));
           if length(betas)~=obj.nCQIs
               error('length of beta calibration parameters must be %d',obj.nCQIs);
           end
           
           betas        = betas(:);        % Store in column format
           obj.betas    = betas;
           obj.betas_dB = 10*log10(betas); % Precalculate also the value in dB
           
           for cqi_ = CQI_range(1):CQI_range(2)
               m_j = CQI_tables(cqi_).modulation_order;
               struct_index = find([BICM_capacity_tables.m_j]==m_j,1,'first');
               BICM_matrix_y(:,cqi_)     = BICM_capacity_tables(struct_index).I;
               BICM_matrix_inv_y(:,cqi_) = interp1(BICM_capacity_tables(struct_index).I_inv,BICM_capacity_tables(struct_index).SNR_inv,BICM_matrix_inv_x);
               
               % Small fix not to have a NaN at the maximum capacity point
               max_capacity = max(BICM_capacity_tables(struct_index).I);
               [C,I] = min(abs(BICM_matrix_inv_x-max_capacity));
               BICM_matrix_inv_y(I,cqi_) = interp1(BICM_matrix_inv_x(1:(I-1)),BICM_matrix_inv_y(1:(I-1),cqi_),BICM_matrix_inv_x(I),'linear','extrap');
           end
           
           obj.length        = length(BICM_matrix_x);
           obj.BICM_matrix_x = BICM_matrix_x;
           obj.BICM_matrix_y = BICM_matrix_y;
           obj.BICM_matrix_inv_x = BICM_matrix_inv_x;
           obj.BICM_matrix_inv_y = BICM_matrix_inv_y;
           obj.multipliers = reshape(0:length(BICM_matrix_x):length(BICM_matrix_x)*(n_CQIs-1),[],1);
           
           % Assume equally-spaced SNRs
           obj.resolution     = (BICM_matrix_x(2)-BICM_matrix_x(1));
           obj.resolution_inv = (BICM_matrix_inv_x(2)-BICM_matrix_inv_x(1));
           
           obj.SINR_range = [min(BICM_matrix_x) max(BICM_matrix_x)];
           obj.BICM_range = [min(BICM_matrix_inv_x) max(BICM_matrix_inv_x)];
           
           % Some (optional) plotting
           if isempty(varargin)
               plot_capacity = false;
           else
               plot_capacity = varargin{1};
           end
           if plot_capacity
               for m_j_idx=1:length(BICM_capacity_tables)
                   BICM_capacity_tables_all(m_j_idx,:) = BICM_capacity_tables(m_j_idx).I;
                   displaynames{m_j_idx} = sprintf('BICM capacity, %d-QAM',2^BICM_capacity_tables(m_j_idx).m_j);
               end
               figure;
               % Assume that all BICM capacities are calcualted over the same
               % SNR range (so I can plot it like this and Matlab automatically
               % puts the colors there)
               plot(BICM_capacity_tables(m_j_idx).SNR,BICM_capacity_tables_all);
               legend(displaynames,'Location','Best');
               xlabel('SNR [dB]');
               ylabel('BICM Capacity I_{m_j}(\gamma)');
               title('BICM capacity');
               grid on;
           end
       end
       
       % Average the give SINR vector. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       % ALL INPUT VALUES IN LINEAR unless an extra parameter is set to "true"
       % Allows for calculations with multiple beta values
       function [effective_SINR_lin effective_SINR_dB] = average(obj,SINR_vector,MCSs,varargin)
           
           if isempty(varargin)
               input_in_dB  = false;
           else
               input_in_dB = varargin{1};
           end
           
           % Put SINR in a row vector
           if size(SINR_vector,1)~=1 && size(SINR_vector,2)==1
               SINR_vector = reshape(SINR_vector,1,[]);
           end
           
           % Data needed beforehand
           SNR_length      = length(obj.BICM_matrix_x);
           multipliers     = obj.multipliers(MCSs);
           multipliers_mat = multipliers(:,ones(length(SINR_vector),1));
           
           if input_in_dB
               SINR_vector_dB = SINR_vector;
           else
               SINR_vector_dB = 10*log10(SINR_vector);
           end
           
           MC_betas_vect_dB    = obj.betas_dB(MCSs);
           SINR_vector_mat_dB  = SINR_vector_dB(ones(length(MCSs),1),:);
           MCS_betas_dB        = MC_betas_vect_dB(:,ones(length(SINR_vector_dB),1));
           SINR_vector_mat_log = SINR_vector_mat_dB - MCS_betas_dB;
           
           SNR_idxs_mat        = round((SINR_vector_mat_log-obj.SINR_range(1))/obj.resolution + 1);
           SNR_idxs_mat(SNR_idxs_mat<1) = 1;
           SNR_idxs_mat(SNR_idxs_mat>SNR_length) = SNR_length;
           
           % Convert to BICM capacity
           Is_mat  = obj.BICM_matrix_y(SNR_idxs_mat+multipliers_mat);
           Is_mean = mean(Is_mat,2);
           
           % Inverse mapping
           Is_idxs = round((Is_mean-obj.BICM_range(1))/obj.resolution_inv + 1);
           Is_idxs(Is_idxs<1) = 1; % Safeguard agains negative indices
           obj_length = obj.length;
           Is_idxs(Is_idxs>obj_length) = obj_length;
           effective_SINR_dB = obj.BICM_matrix_inv_y(Is_idxs+multipliers)+MC_betas_vect_dB; % Version with scaling
           effective_SINR_lin = 10.^(effective_SINR_dB/10);
       end
   end
end 
