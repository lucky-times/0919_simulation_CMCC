function activate_UE_id = NR_get_activate_UE( UEs )
% 获得需要进行统计的UE的id
% 输出参数
% activate_UE_id：需要进行统计的UE的id

activate_UE_id = [];
for u_ = 1:length(UEs)
    if ~UEs(u_).deactivate_UE
        if isempty(activate_UE_id)
            activate_UE_id = UEs(u_).id;
        else
            activate_UE_id = [activate_UE_id UEs(u_).id];
        end
    end
end

end

