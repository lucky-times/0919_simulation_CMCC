function [neighboring_eNodeBs neighboring_eNodeBs_sectors] = LTE_init_get_eNodeB_neighbors(this_eNodeB,eNodeBs,varargin)
% Generates a list of neighboring eNodeBs for the given eNodeB.
% If distance<=distance_threshold, then 2 given eNodeBs are neighbors.
% (c) Josep Colom Ikuno, Martin Wrulich INTHFT, 2008
% input:    this_eNodeB          ... our considered eNodeB
%           eNodeBs              ... all the eNodeBs. Note that
%                                    eNodeBs(i).id==i must be fulfilled!!
% output:   neighboring_eNodeBs  ... the neighboring eNodeBs

% We will store the distances here
% eNodeB_ID, distance

if ~isempty(varargin)
    % Take all sites closer than distance_threshold
    distance_mode = true;
    distance_threshold = varargin{1};
else
    % Take the 6 closest sites
    distance_mode = false;
end

if length(eNodeBs)==1
    neighboring_eNodeBs = [];
else
    allSites = reshape([eNodeBs.pos],2,[])';
    allDistances = sqrt(sum((repmat(this_eNodeB.pos,[size(allSites,1) 1])-allSites).^2,2));
    if distance_mode
        neighbors_idx = allDistances<=distance_threshold&allDistances~=0;
        neighboring_eNodeBs = eNodeBs(neighbors_idx);
    else
        [allDistances_sorted indexes] = sortrows(allDistances);
        to_take = 6;
        lowest_idxs = indexes(2:min(length(indexes),2+to_take-1));
        lowest_idxs = sortrows(lowest_idxs);
        neighboring_eNodeBs = eNodeBs(lowest_idxs);
    end
end

if ~isempty(neighboring_eNodeBs)
    neighboring_eNodeBs_sectors_idx = 1; % Initialization
    for b_=1:length(neighboring_eNodeBs)
        for s_=1:length(neighboring_eNodeBs(b_).sectors)
            if neighboring_eNodeBs_sectors_idx==1
                neighboring_eNodeBs_sectors = neighboring_eNodeBs(b_).sectors(s_);
            else
                neighboring_eNodeBs_sectors(neighboring_eNodeBs_sectors_idx) = neighboring_eNodeBs(b_).sectors(s_);
            end
            neighboring_eNodeBs_sectors_idx = neighboring_eNodeBs_sectors_idx + 1;
        end
    end
else
    neighboring_eNodeBs_sectors = [];
end

