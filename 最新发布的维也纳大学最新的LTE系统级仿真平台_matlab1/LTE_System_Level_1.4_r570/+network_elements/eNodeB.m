classdef eNodeB < handle
    % Class that represents an LTE eNodeB
    % Attributes of this give eNodeB
    % (c) Josep Colom Ikuno, INTHFT, 2008
    properties
        % eNodeB ID. It is also the index in the eNodeBs array in the
        % simulator!!!
        id
        % pos in meters (x,y)
        pos
        % stores info about the eNodeB's sectors. Mainly azimuth and antenna type
        sectors
        % the neighboring sites, ordered by distance
        neighbors
        % the neighboring eNodeBs (i.e. the eNodeBs in the nighboring sites), ordered by distance
        neighbors_eNodeB
        % Whether this enodeB is the target cell
        target_cell = false;
        
        % extra params that may not be always used (mainly extra information from the network planning tool)
        name
        altitude = 0;
        site_name
        
        % network clock. Tells the eNodeB in what TTI he is
        clock
    end
    % Associated eNodeB methods
    methods
        function print(obj)
            fprintf('eNodeB %d, position (x,y): (%3.2f,%3.2f), altitude: %3.0fm, %d attached UEs\n',obj.id,obj.pos(1),obj.pos(2),obj.altitude,obj.attached_UEs);
            fprintf('  Neighbor eNodeBs: ');
            for n_=1:length(obj.neighbors)
                fprintf('%d ',obj.neighbors(n_).id);
            end
            fprintf('\n');
            for s_=1:length(obj.sectors)
                obj.sectors(s_).print;
            end
        end
        % Attachs a user to this eNodeB and sector, first checking that the node is
        % not already in the list
        function attachUser(obj,user,sector)
            obj.sectors(sector).attachUser(user);
        end
        % Deattaches a user from this eNodeB.
        function deattachUser(obj,user)
            obj.sectors(user.attached_sector).deattachUser(user);
        end
        % Queries whether a user is attached
        function is_attached = userIsAttached(obj,user)
            for s_ = 1:length(obj.sectors)
                if obj.sectors(1).userIsAttached(user);
                    is_attached = true;
                    return
                end
            end
            is_attached = false;
        end
        % Returns the number of UEs currently attached to this eNodeB
        function number_of_atached_UEs = attached_UEs(obj)
            temp = zeros(1,length(obj.sectors));
            for s_ = 1:length(obj.sectors)
                temp(s_) = obj.sectors(s_).attached_UEs;
            end
            number_of_atached_UEs = sum(temp);
        end
    end
end