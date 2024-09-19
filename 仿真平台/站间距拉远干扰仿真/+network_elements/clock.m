classdef clock < handle
% 时钟类，所有network element 共享一个时钟

   properties
       current_TTI  % TTI数
       TTI_time     % 一个TTI持续的时间
       time         % 实际时间
   end

   methods
       % 构造函数
       function obj = clock(TTI_time)
           obj.current_TTI = 0;
           obj.time        = 0;
           obj.TTI_time    = TTI_time;
       end
       
       % TTI加一
       function advance_1_TTI(obj)
           obj.current_TTI = obj.current_TTI + 1;
           obj.time = obj.time + obj.TTI_time;
       end
   end
end 
