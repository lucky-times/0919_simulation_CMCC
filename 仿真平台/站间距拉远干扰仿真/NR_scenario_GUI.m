function varargout = NR_scenario_GUI(varargin)
% NR_SCENARIO_GUI MATLAB code for NR_scenario_GUI.fig
%      NR_SCENARIO_GUI, by itself, creates a new NR_SCENARIO_GUI or raises the existing
%      singleton*.
%
%      H = NR_SCENARIO_GUI returns the handle to a new NR_SCENARIO_GUI or the handle to
%      the existing singleton*.
%
%      NR_SCENARIO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NR_SCENARIO_GUI.M with the given input arguments.
%
%      NR_SCENARIO_GUI('Property','Value',...) creates a new NR_SCENARIO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NR_scenario_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NR_scenario_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NR_scenario_GUI

% Last Modified by GUIDE v2.5 23-Jun-2017 13:51:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NR_scenario_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NR_scenario_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NR_scenario_GUI is made visible.
function NR_scenario_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NR_scenario_GUI (see VARARGIN)

% Choose default command line output for NR_scenario_GUI
% 界面初始化函数
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.homo,'Value',1); % 初始化默认为同构情景
set(handles.UMA,'Value',1); % 初始化默认为UMA场景
set(handles.ue_bf,'Value',1); % 初始化默认UE采用BF
set(handles.bs_bf,'Value',1); % 默认BS使用BF
set(handles.issave,'Value',1); % 默认保存波束赋型数据
set(handles.use_cache,'Value',1); % 默认使用波束赋型缓存
set(handles.still,'Value',1); % 默认UE是静止的
set(handles.wrap,'Value',1); %默认不使用wraparound

set(handles.Rdrop,'Enable','off'); % 默认同构系统为UMA系统，不允许修改UMI的类型
set(handles.femto,'Enable','off'); % 默认同构系统为UMA系统，不允许修改InH的类型
set(handles.UMA2UMI,'Enable','off'); % 不允许选择UMA_to_UMI的异构场景
set(handles.UMI2InH,'Enable','off'); % 不允许选择UMI_to_InH的异构场景
set(handles.UMA2InH,'Enable','off'); % 不允许选择UMA_to_InH的异构场景
set(handles.F_2_L,'Enable','off'); % 不允许选择前者干扰后者
set(handles.L_2_F,'Enable','off'); % 不允许选择后者干扰前者
set(handles.bs_max_gain2,'Enable','off'); % 不允许填入异构最大单元增益
set(handles.h_bs2,'Enable','off'); % 不允许填入异构基站高度
set(handles.p_bs2,'Enable','off'); % 不允许填入异构基站发射功率
set(handles.speed,'Enable','off'); % 不允许填入速度（因为用户移动没有选择）
set(handles.TTI,'Enable','off'); % 不允许修改仿真时长（因为用户移动没有选择）
set(handles.div_angle,'Enable','off'); % 不允许填入偏移角度
set(handles.div_dist,'Enable','off'); % 不允许填入偏移距离


set(handles.bs_array_out,'String','[8,16]'); % 默认户外天线阵配置为8*16
set(handles.bs_array_in,'String','[4,8]'); % 默认室内天线阵配置为4*8
set(handles.ue_max_gain,'String','5'); % 默认用户天线最大单元增益为5dB
set(handles.bs_max_gain1,'String','8'); % 默认victim系统（和同构下aggressor系统）基站天线最大单元增益为8dB
set(handles.bs_phi,'String','0'); % 默认基站侧波束赋水平角误差上界为0
set(handles.bs_theta,'String','0'); % 默认基站侧波束垂直角误差上界为0
set(handles.ue_phi,'String','0'); % 默认UE侧波束水平角误差上界为0
set(handles.ue_theta,'String','0'); % 默认UE侧波束垂直角误差上界为0
set(handles.acir_lower,'String','5'); % 默认ACIR研究下界为5dB
set(handles.acir_upper,'String','40'); % 默认ACIR研究上界为40dB
set(handles.acir_step,'String','5'); % 默认ACIR研究步进值为5dB
set(handles.acir,'String','15'); % 默认ACIR校准值为15dB
set(handles.slot,'String','5'); % 默认快照数为10
set(handles.f1,'String','30'); % 默认victim系统频率为30GHz
set(handles.f2,'String','30.2'); % 默认aggressor系统频率为30.2GHz
set(handles.band,'String','200'); % 默认带宽为200MHz
set(handles.p_ue,'String','23'); % 默认UE的发射功率为23dBm
set(handles.p_bs1,'String','43'); % 默认victim系统（和同构下aggressor系统）基站发射功率为43dBm
set(handles.nf_ue,'String','11'); % 默认UE侧噪声系数为11dB
set(handles.nf_bs,'String','11'); % 默认BS侧噪声系数为11dB
set(handles.h_bs1,'String','25'); % 默认victim系统（和同构下aggressor系统）基站高度为25m
set(handles.shadow_r,'String','0.5'); % 默认阴影衰落相关系数为0.5
set(handles.res,'String','5'); % 默认分辨率为5
set(handles.isd,'String','200'); % 默认系统ISD为200m
set(handles.speed,'String','0'); %默认用户速度为0
set(handles.TTI,'String','1'); % 默认仿真时间为1个TTI
set(handles.UEs_per_slot,'String','1'); % 默认每次快照服务1个用户


% UIWAIT makes NR_scenario_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NR_scenario_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in structure.
function structure_Callback(hObject, eventdata, handles)
% hObject    handle to structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns structure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from structure


% --- Executes during object creation, after setting all properties.
function structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hete.
function hete_Callback(hObject, eventdata, handles)
% hObject    handle to hete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% 异构按钮触发器

set(handles.homo,'Value',0); %空置同构选项
%空置具体同构场景的选项
set(handles.RMA,'Value',0);  
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
set(handles.UMI,'Value',0);
set(handles.F_2_L,'Enable','on');%允许选择正向干扰
set(handles.L_2_F,'Enable','on');%允许选择反向干扰
set(handles.F_2_L,'Value',1);%默认异构场景是前者干扰后者（正向）
set(handles.L_2_F,'Value',0);
%默认异构场景为UMA_to_UMI
set(handles.UMA2UMI,'Value',1);
set(handles.loc,'Value',2); %异构场景默认选择不共址
%锁死具体同构场景的选项
set(handles.RMA,'Enable','off');
set(handles.UMA,'Enable','off');
set(handles.InH,'Enable','off');
set(handles.UMI,'Enable','off');
set(handles.Rdrop,'Enable','on'); % 默认异构系统为UMa_to_UMi.故解锁UMI特殊场景选择，支持UMA与manhattan和UMA与Rdrop两种类型
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','off'); % 锁死InH类型选项，在此场景下无用
set(handles.wrap,'Enable','off'); %锁死wraparound选项
% 具体异构场景解锁
set(handles.UMA2UMI,'Enable','on'); 
set(handles.UMA2InH,'Enable','on');
set(handles.UMI2InH,'Enable','on');
set(handles.bs_max_gain2,'Enable','on'); %解锁异构系统基站天线最大单元增益的设置
set(handles.h_bs2,'Enable','on'); %解锁异构系统基站高度的设置 
set(handles.p_bs2,'Enable','on'); %解锁异构系统基站发射功率的设置
set(handles.loc,'Enable','off'); %锁死共址的选择，禁止选择为共址选项
% set(handles.together,'Enable','off'); %锁死同步/异步干扰类型——异构系统不支持异步干扰的情况
set(handles.div_angle,'Enable','off'); %解锁小场景偏离大场景中心的角度
set(handles.div_dist,'Enable','off'); %解锁小场景偏离大场景中心的距离
%默认小场景在大场景的中心
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.sysnum,'Enable','off'); %锁死系统数选项——异构系统的系统数为2
set(handles.isd,'Enable','on'); %解锁ISD设置，允许自行设置victim系统的ISD
%解锁ACIR相关的设置——单系统选项时会将ACIR相关值清空并且锁死
set(handles.acir_upper,'Enable','on'); 
set(handles.acir_lower,'Enable','on');
set(handles.acir_step,'Enable','on');
set(handles.acir,'Enable','on');
%锁死系统类型选择——只有UMA/RMA支持LTE与NR的共存干扰
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','on');


set(handles.bs_max_gain2,'String','8'); %触发设置异构基站天线最大单元增益为8dB
set(handles.h_bs2,'String','25'); %触发设置异构系统基站高度为25m
set(handles.p_bs2,'String','43'); %触发设置异构系统基站发射功率为43dBm
set(handles.p_bs1,'String','33'); %触发设置victim系统基站发射功率为33dBm
set(handles.h_bs1,'String','10'); %触发设置victim系统基站高度为10m
set(handles.slot,'String','5'); %触发设置系统快照数为5
 %重置两系统频率为30G、30.2GHz——RMa场景会触发频率的更改
set(handles.f1,'String','30'); 
set(handles.f2,'String','30.2');
%触发设置系统采用同步干扰模式
set(handles.together,'Value',1);
set(handles.sysnum,'Value',1); %触发设置系统数为2
set(handles.wrap,'Value',1); %触发设置不采用wraparound——只有UMA/RMA场景支持wraparound
set(handles.isd,'String','200'); %触发重置系统ISD我200
%触发设置ACIR相关参数的设置——如果从单系统状态触发到异构状态，则需要解锁ACIR相关参数并赋值
set(handles.acir,'String','15');
set(handles.acir_step,'String','5');
set(handles.acir_lower,'String','5');
set(handles.acir_upper,'String','40');
set(handles.handover,'Value',1); % 触发设置采用3dB handover
set(handles.shadow_r,'String','0.5'); % 触发设置阴影衰落相关系数为0.5
%重置两系统的类型为NR
set(handles.a_type,'Value',1);
set(handles.v_type,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of hete


% --- Executes on button press in homo.
function homo_Callback(hObject, eventdata, handles)
% hObject    handle to homo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 点击同构选框的触发事件
% 默认为UMA系统，故锁死两种场景的类型
set(handles.Rdrop,'Enable','off');
set(handles.femto,'Enable','off');
%锁死异构系统具体场景的选择并空置值
set(handles.UMA2UMI,'Value',0);
set(handles.UMI2InH,'Value',0);
set(handles.UMA2InH,'Value',0);
set(handles.UMA2UMI,'Enable','off');
set(handles.UMI2InH,'Enable','off');
set(handles.UMA2InH,'Enable','off');
set(handles.F_2_L,'Enable','off');%不允许选择正向干扰
set(handles.L_2_F,'Enable','off');%不允许选择反向干扰
%异构干扰方向均复位
set(handles.F_2_L,'Value',0);
set(handles.L_2_F,'Value',0);
%打开同构系统具体场景的选择
set(handles.UMA,'Enable','on');
set(handles.RMA,'Enable','on');
set(handles.UMI,'Enable','on');
set(handles.InH,'Enable','on');
%锁死异构系统专有参数（最大单元增益，天线高度，发射功率）
set(handles.bs_max_gain2,'Enable','off');
set(handles.h_bs2,'Enable','off');
set(handles.p_bs2,'Enable','off');
%解锁异步干扰选择
set(handles.together,'Enable','on');
set(handles.sysnum,'Enable','on'); %解锁单系统选项
%由于异构系统默认选择UMA场景
set(handles.loc,'Enable','on'); %解锁不共址选项
set(handles.wrap,'Enable','on'); %解锁wraparound的选项


set(handles.bs_array_out,'String','[8,16]');
set(handles.ue_max_gain,'String','5');
set(handles.bs_max_gain1,'String','8'); %重置基站侧天线最大单元增益
set(handles.slot,'String','5');% 重置快照数
set(handles.h_bs1,'String','25'); % 重置基站侧天线高度
set(handles.res,'String','5'); %重置分辨率
set(handles.isd,'String','200'); %重置ISD
%清空异构系统专有参数
set(handles.h_bs2,'String','');  
set(handles.p_bs2,'String','');
set(handles.bs_max_gain2,'String','');

set(handles.together,'Value',1); %默认为同步干扰
set(handles.hete,'Value',0); %空置同构选项——同构异构只能选择一个
set(handles.sysnum,'Value',1); %重置系统数为双系统
set(handles.UMA,'Value',1); %默认同构场景为UMA
set(handles.loc,'Value',1); %重置共址选项
%清空两系统偏置角度与偏置值
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.wrap,'Value',1); %默认不采用wraparound
%解锁偏置角度与距离——默认同构系统为UMA
set(handles.v_type,'Enable','on'); 
set(handles.a_type,'Enable','on');

% Hint: get(hObject,'Value') returns toggle state of homo


% --- Executes on selection change in sysnum.
function sysnum_Callback(hObject, eventdata, handles)
% hObject    handle to sysnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%操作单/双系统操作栏

n_sysnum = get(hObject,'value'); %得到选择的内容，1为双系统，2为单系统
if n_sysnum==1 %双系统
    set(handles.loc,'Value',1); %默认共址情景
    if get(handles.UMA,'Value') + get(handles.RMA,'Value')>0 %如果采用UMA或RMA
        set(handles.loc,'Enable','on'); %解锁不共址选项 ——只有UMA或RMA支持自行设置偏置状态
    else
        set(handles.loc,'Enable','off'); %锁死不共址选项
    end
    %解锁ACIR相关设置并恢复ACIR的相关默认设置
    set(handles.acir_lower,'String','5');
    set(handles.acir_upper,'String','40');
    set(handles.acir,'String','15');
    set(handles.acir_step,'String','5');
    set(handles.acir_upper,'Enable','on');
    set(handles.acir_lower,'Enable','on');
    set(handles.acir_step,'Enable','on');
    set(handles.acir,'Enable','on');
    %解锁干扰模式选择并默认为同步干扰
    set(handles.together,'Value',1);
    set(handles.together,'Enable','on');
else %单系统
    %锁死并空置共址选项
    set(handles.loc,'Value',0);
    set(handles.loc,'Enable','off');
    %锁死并空置偏置角度与偏置距离
    set(handles.div_angle,'Enable','off');
    set(handles.div_dist,'Enable','off');
    set(handles.div_angle,'String','');
    set(handles.div_dist,'String','');
    %锁死并空置ACIR的相关设置
    set(handles.acir_lower,'String','');
    set(handles.acir_upper,'String','');
    set(handles.acir,'String','');
    set(handles.acir_step,'String','');
    set(handles.acir_upper,'Enable','off');
    set(handles.acir_lower,'Enable','off');
    set(handles.acir_step,'Enable','off');
    set(handles.acir,'Enable','off');
    %锁死干扰模式的选择，并重置为同步干扰
    set(handles.together,'Enable','off');
    set(handles.together,'Value',0);
end
   
% Hints: contents = cellstr(get(hObject,'String')) returns sysnum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sysnum


% --- Executes during object creation, after setting all properties.
function sysnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sysnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)

% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 点击确定按钮的触发事件
% 由于不能想z_test那样清空所有的变量（会清空GUI的句柄），只能采用clc清空部分数据。如图形界面运行出错而脚本运行正常，
% 则运行clear后再运行GUI
% 该触发事件主要是进行参数的传递和主函数，拓扑GUI的调用
% 从该函数中可以很清楚地看出句柄变量与参数的对应关系
clc;
SYS_config = NR_load_params;%导入默认仿真参数
%% 灵活性较弱的变量，默认设置，不在图形界面中配置
SYS_config.shadow_fading_type         = 'claussen';
SYS_config.compact_results_file       = true;
SYS_config.TTI_length = 1;% Length of a TTI (subframe), in seconds.
SYS_config.default_shown_GUI_cells = [];
SYS_config.PCe_param = 0; %功控的精确度
SYS_config.interference_type = 2;% 0同频干扰，1邻频干扰，2都有，只在双系统下有效
SYS_config.attatch_mode=3;  % 1表示路损加BS端单元增益；2表示路损加UE端单元增益；3表示路损加BS和UE端单元增益。用来确定UE归属小区（撒点）
SYS_config.Gama = 1; %功控中的gamma参数
SYS_config.antenna.antenna_gain_pattern = 'NRAntennaBeamforming'; %波束赋型的类型（可以在此处增加扩展其他的波束赋型图案）
SYS_config.AntennaPattern3d = true;
SYS_config.beam_loss=0; %波束损耗
SYS_config.cable_loss=0; %天线馈线损耗

SYS_config.ISD = str2double(get(handles.isd,'String'));
%% 专用参数
if get(handles.UMA,'Value')==1 %UMA
    SYS_config.scene_type = 'UMA';
    SYS_config.isNewUMA = true;
    SYS_config.sector_azimuths = 0:360/3:359;
    SYS_config.shadow_fading_sd_LOS = 4; %LOS下阴影衰落标准差
    SYS_config.shadow_fading_sd_NLOS = 6; %NLOS下阴影衰落标准差
    SYS_config.default_shown_GUI_cells = 1:57;
    SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro'; %第一个系统的环境
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro'; %第二个系统的环境
end
if get(handles.UMI,'Value')==1 %UMI
    SYS_config.scene_type = 'UMI';
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths = 0;% 初始化扇区角度
    SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_NLOS = 7.82;
end
if get(handles.RMA,'Value')==1 %RMA
    SYS_config.scene_type = 'RMa';
    SYS_config.sector_azimuths = 0:360/3:359;
    SYS_config.macroscopic_pathloss_model_settings.environment = 'rural_macro';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'rural_macro';
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_LOS2 = 6; %900中RMA场景下LOS状态有两个LOS。
    SYS_config.shadow_fading_sd_NLOS = 8;
    SYS_config.default_shown_GUI_cells = 1:57;
end
if get(handles.InH,'Value')==1 %InH
    SYS_config.scene_type = 'InH';
    SYS_config.sector_azimuths = 0;% 初始化扇区角度
    SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.UE_height=1;
end
if get(handles.UMA2UMI,'Value')==1 %UMa_to_UMi
    SYS_config.scene_type = 'UMa_to_UMi';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.isNewUMA = true;
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths = 0;% 初始化扇区角度
    if SYS_config.isS2F %如果是前者干扰后者，则为干扰系统配置UMA场景，被干扰系统配置INH场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
    else %否则为干扰系统配置indoor场景，为被干扰系统配置UMA场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    end
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_NLOS = 7.82;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 6;
    SYS_config.sector_azimuths2 = 0:360/3:359;
   
end
if get(handles.UMA2InH,'Value')==1 %UMa_to_InH
    SYS_config.scene_type = 'UMa_to_InH';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.sector_azimuths = 0;% 初始化扇区角度
    if SYS_config.isS2F %如果是前者干扰后者，则为干扰系统配置UMA场景，被干扰系统配置INH场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
    else %否则为干扰系统配置indoor场景，为被干扰系统配置UMA场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    end
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.sector_azimuths2 = 0:360/3:359;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 6;
    SYS_config.UE_height=1.5;
end
if get(handles.UMI2InH,'Value')==1 %UMi_to_InH
    SYS_config.scene_type = 'UMi_to_InH';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.sector_azimuths = 0;% 初始化扇区角度
    if SYS_config.isS2F %如果是前者干扰后者，则为干扰系统配置UMI场景，被干扰系统配置INH场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    else %否则为干扰系统配置indoor场景，为被干扰系统配置UMI场景
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    end
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths2 = 0;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 7.82;
    SYS_config.UE_height=1.5;
end
if get(handles.Rdrop,'Value')==1
    SYS_config.isManhattan = false;
else
    SYS_config.isManhattan = true;
end
if get(handles.femto,'Value')==1
    SYS_config.isFemto = false;
else
    SYS_config.isFemto = true;
end

SYS_config.simulation_time_tti = str2double(get(handles.TTI,'String'));
SYS_config.keep_UEs_still  = get(handles.still,'Value')==1;
SYS_config.UE_speed = str2double(get(handles.TTI,'String'))/3.6; %将速度值转换为m/s

SYS_config.beam_accuracy_phi_ue = str2double(get(handles.ue_phi,'String'));
SYS_config.beam_accuracy_theta_ue = str2double(get(handles.ue_theta,'String'));
SYS_config.beam_accuracy_phi_bs = str2double(get(handles.bs_phi,'String'));
SYS_config.beam_accuracy_theta_bs = str2double(get(handles.bs_theta,'String'));

if get(handles.distribution,'Value')==1
    SYS_config.beam_accurancy_type = 'Constant';
else
    SYS_config.beam_accurancy_type = 'Unique';
end

SYS_config.asynchronization = get(handles.together,'Value')==2;


SYS_config.frequency = str2double(get(handles.f1,'String')) * 1e9; %单位转换，从GHz转换为Hz
SYS_config.frequency2 = str2double(get(handles.f2,'String')) * 1e9; %单位黄钻换，从GHz转换为Hz
SYS_config.bandwidth = str2double(get(handles.band,'String')) * 1e6; %单位转换，从MHz转换为Hz
SYS_config.photo = get(handles.isphoto,'Value');
SYS_config.ACIR_dB = str2double(get(handles.acir,'String'));
SYS_config.isACIR_loop = true;
SYS_config.ACIR_lower = str2double(get(handles.acir_lower,'String'));
SYS_config.ACIR_upper = str2double(get(handles.acir_upper,'String'));
SYS_config.isDouble = get(handles.sysnum,'Value')==1;
if get(handles.ue_bf,'value')==1 && get(handles.bs_bf,'value')==1
    SYS_config.antenna_mode = 2;
elseif get(handles.bs_bf,'value')==1
    SYS_config.antenna_mode = 1;
else
    SYS_config.antenna_mode = 0;
end
SYS_config.use_cache = get(handles.use_cache,'value')==1;
SYS_config.isSave = get(handles.issave,'value')==1;
SYS_config.antenna_element = str2num(get(handles.bs_array_out,'String'));
SYS_config.antenna_element_InH = str2num(get(handles.bs_array_in,'String'));
SYS_config.UE_tx_power = 10^(0.1*str2double(get(handles.p_ue,'String')))/1000;
SYS_config.UE_max_antenna_gain = str2double(get(handles.ue_max_gain,'String'));
SYS_config.UE_receiver_noise_figure = str2double(get(handles.nf_ue,'String'));
SYS_config.BS_receiver_noise_figure = str2double(get(handles.nf_bs,'String'));
if get(handles.v_type,'Value')==1
    SYS_config.macroscopic_pathloss_model = 'TS38901'; %901对应NR系统的路损模型
else
    SYS_config.macroscopic_pathloss_model = 'TS36942'; %942对应LTE系统的路损模型
end
if get(handles.a_type,'Value')==1
    SYS_config.macroscopic_pathloss_model2 = 'TS38901';
else
    SYS_config.macroscopic_pathloss_model2 = 'TS36942';
end
SYS_config.isWraparound = get(handles.wrap,'Value')==2;
SYS_config.ISD = str2double(get(handles.isd,'String'));
if get(handles.loc,'Value')==1
    SYS_config.shift_mode = 0;
else
    SYS_config.shift_mode = 2;
%     if get(handles.UMA,'Value') + get(handles.RMA,'Value')==0 
%         SYS_config.shift_mode = 1;
%     else
%         SYS_config.shift_mode = 2;
%     end
end
SYS_config.angle_between_systems = str2double(get(handles.div_angle,'String'));
SYS_config.ISD_between_systems = str2double(get(handles.div_dist,'String'));
SYS_config.hand_over = get(handles.handover,'Value')==1;
SYS_config.map_resolution = str2double(get(handles.res,'String'));
SYS_config.shadow_fading_map_resolution = str2double(get(handles.res,'String'));
SYS_config.n_snapshot = str2double(get(handles.slot,'String'));
SYS_config.n_snapshot2 = str2double(get(handles.slot,'String'));
SYS_config.n_UE_served_per_BS = str2double(get(handles.UEs_per_slot,'String'));
SYS_config.r_eNodeBs = str2double(get(handles.shadow_r,'String'));
SYS_config.site_height = str2double(get(handles.h_bs1,'String'));
SYS_config.eNodeB_tx_power = 10^(0.1*str2double(get(handles.p_bs1,'String')))/1000; %将dBm转换为W
SYS_config.antenna.max_antenna_gain = str2double(get(handles.bs_max_gain1,'String')); %将字符转换为数值
SYS_config.eNodeB_tx_power2 = 10^(0.1*str2double(get(handles.p_bs2,'String')))/1000; %将dBm转换为W
SYS_config.site_height2 = str2double(get(handles.h_bs2,'String'));
SYS_config.antenna.max_antenna_gain2 = str2double(get(handles.bs_max_gain2,'String'));

 
 output_results_file = NR_sim_main(SYS_config); %运行主程序
 simulation_data                   = load(output_results_file); %导出数据
 
 GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data); %绘制拓扑图
% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 单击取消按钮触发的时间
clear; %清空所有的变量
close(gcf); %关闭GUI

% --- Executes on button press in UMA2InH.
function UMA2InH_Callback(hObject, eventdata, handles)
% hObject    handle to UMA2InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 触发UMA干扰InH场景的配置
% InH类型可选，不涉及UMI场景，故其类型不可选
set(handles.Rdrop,'Enable','off');
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','on');
% 空置其他异构场景的选中状态
set(handles.UMA2UMI,'Value',0);
set(handles.UMI2InH,'Value',0);
%重置两系统基站参数
set(handles.p_bs1,'String',23);
set(handles.h_bs1,'String',3);
set(handles.p_bs2,'String',43);
set(handles.h_bs2,'String',25);
set(handles.bs_max_gain1,'String',5);
set(handles.res,'String',2); %重置分辨率
set(handles.slot,'String',5); %重置快照数
%重置阴影衰落相关系数
set(handles.shadow_r,'String','0.5');
set(handles.v_type,'Value',1);%默认被干扰系统为NR
set(handles.a_type,'Value',1);%默认干扰系统为NR
if get(handles.F_2_L,'Value') %如果干扰系统是UMA
    set(handles.v_type,'Enable','off'); %被干扰系统一定为NR且不可变
    set(handles.a_type,'Enable','on'); %干扰系统可选LTE
else  %如果被干扰系统是UMA
    set(handles.v_type,'Enable','on');%被干扰系统可选LTE
    set(handles.a_type,'Enable','off');%干扰系统一定为NR
end
%解锁偏置方向角和偏置距离
set(handles.div_angle,'Enable','on');
set(handles.div_dist,'Enable','on');
set(handles.div_angle,'String','0');
set(handles.div_dist,'String','0');
% Hint: get(hObject,'Value') returns toggle state of UMA2InH


% --- Executes on button press in UMI2InH.
function UMI2InH_Callback(hObject, eventdata, handles)
% hObject    handle to UMI2InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 触发UMI干扰InH场景的配置
% 由于不支持manhattan与InH的异构场景，故该场景下InH类型可选，UMI场景不可选
set(handles.Rdrop,'Value',1); %重置为random drop场景研究与InH场景的共存
set(handles.Rdrop,'Enable','off');
set(handles.femto,'Enable','on');
% 空置其他异构场景的选中状态
set(handles.UMA2UMI,'Value',0);
set(handles.UMA2InH,'Value',0);
%重置两系统基站参数
set(handles.p_bs1,'String',23);
set(handles.h_bs1,'String',3);
set(handles.p_bs2,'String',33);
set(handles.h_bs2,'String',10);
set(handles.bs_max_gain1,'String',5);
set(handles.res,'String',2); %重置分辨率
set(handles.slot,'String',5); %重置快照数
set(handles.shadow_r,'String','0.5'); %重置阴影衰落相关系数
%两系统均为NR且不可选LTE
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of UMI2InH
%解锁偏置方向角和偏置距离
set(handles.div_angle,'Enable','on');
set(handles.div_dist,'Enable','on');
set(handles.div_angle,'String','0');
set(handles.div_dist,'String','0');

% --- Executes on button press in UMA2UMI.
function UMA2UMI_Callback(hObject, eventdata, handles)
% hObject    handle to UMA2UMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 触发UMA干扰UMI场景的配置、
% 平台支持UMA与manhattan、UMA与Rdrop共存的两种异构场景，故UMI类型可选，不涉及InH场景，故无需选择InH场景
set(handles.Rdrop,'Enable','on');
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','off');
% 空置其他异构场景的选中状态
set(handles.UMA2InH,'Value',0);
set(handles.UMI2InH,'Value',0);
%重置两系统基站参数
set(handles.p_bs1,'String',33);
set(handles.h_bs1,'String',10);
set(handles.p_bs2,'String',43);
set(handles.h_bs2,'String',25);
set(handles.bs_max_gain1,'String',8);
set(handles.res,'String',5); %重置分辨率
set(handles.slot,'String',5);%重置快照数
 %重置阴影衰落相关系数
set(handles.shadow_r,'String','0.5');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
if get(handles.F_2_L,'Value') %如果干扰系统是UMA
    set(handles.v_type,'Enable','off'); %被干扰系统一定为NR且不可变
    set(handles.a_type,'Enable','on'); %干扰系统可选LTE
else  %如果被干扰系统是UMA
    set(handles.v_type,'Enable','on');%被干扰系统可选LTE
    set(handles.a_type,'Enable','off');%干扰系统一定为NR
end
%锁死偏置方向角和偏置距离
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
    
% Hint: get(hObject,'Value') returns toggle state of UMA2UMI


% --- Executes on button press in InH.
function InH_Callback(hObject, eventdata, handles)
% hObject    handle to InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 单击选中InH场景下的触发事件
% 空置其他同构场景的选中状态
set(handles.RMA,'Value',0);
set(handles.UMA,'Value',0);
set(handles.UMI,'Value',0);
%锁死wraparound选项并重置为不采用wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','off');
set(handles.Rdrop,'Enable','off');%锁死UMI场景的类型
set(handles.femto,'Enable','on');%解锁InH场景的类型
set(handles.slot,'String','5'); %重置快照数
%重置基站相关参数
set(handles.p_bs1,'String','23');
set(handles.h_bs1,'String','3');
set(handles.bs_max_gain1,'String','5');
%清空并锁死偏置状态的设置
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','0.5'); %重置分辨率
%重置频率的设置
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%锁死并重置ISD设置——InH场景不允许修改ISD
set(handles.isd,'String','20');
set(handles.isd,'Enable','off');
set(handles.handover,'Value',2); %重置handover状态为不采用3dB handover
set(handles.shadow_r,'String','0'); %重置阴影衰落的相关系数为0
if get(handles.sysnum,'Value')==1 %如果是双系统则关闭共址的选择，设置默认值为共址——InH不支持系统偏置
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','off');
end
% set(handles.sysnum,'Value',1);
% 锁死并重置系统类型选择
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);

% Hint: get(hObject,'Value') returns toggle state of InH


% --- Executes on button press in UMI.
function UMI_Callback(hObject, eventdata, handles)
% hObject    handle to UMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 单击选中UMi场景下的触发事件
% 空置其他同构场景的选中状态
set(handles.RMA,'Value',0);
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
%锁死wraparound选项并重置为不采用wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','off');
%解锁UMI场景类型并默认设为Random drop类型
set(handles.Rdrop,'Enable','on');
set(handles.femto,'Enable','off');% 锁死InH类型选择
set(handles.Rdrop,'Value',1);
%重置基站相关参数
set(handles.p_bs1,'String','33');
set(handles.slot,'String','5');
set(handles.h_bs1,'String','10');
set(handles.bs_max_gain1,'String','8');
%清空并锁死偏置状态的设置
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','5'); %重置分辨率
%重置频率的设置
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%锁死并重置ISD设置——InH场景不允许修改ISD
set(handles.isd,'String','200');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1); %重置handover状态为采用3dB handover
if get(handles.sysnum,'Value')==1 %如果是双系统则关闭共址的选择，设置默认值为不共址——UMI为不共址情形
    set(handles.loc,'Value',2);
    set(handles.loc,'Enable','off');
end
set(handles.shadow_r,'String','0.5');%重置阴影衰落的相关系数为0.5
% 锁死并重置系统类型选择
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
% set(handles.sysnum,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of UMI


% --- Executes on selection change in loc.
function loc_Callback(hObject, eventdata, handles)
% hObject    handle to loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 当选择共址状态栏时触发的事件 1为共址 2为不共址（只有同构时可以修改）
loc_type = get(hObject,'value');
if loc_type==2 %不共址
    if get(handles.UMA,'Value') + get(handles.RMA,'Value')>0 %RMA或者UMA则解锁并设置偏置状态的默认值
        set(handles.div_angle,'Enable','on');
        set(handles.div_dist,'Enable','on');
        set(handles.div_angle,'String','45');
        set(handles.div_dist,'String','50');
    end
else %共址则锁死并清空偏置状态
    set(handles.div_angle,'Enable','off');
    set(handles.div_dist,'Enable','off');
    set(handles.div_angle,'String','');
    set(handles.div_dist,'String','');
end

% Hints: contents = cellstr(get(hObject,'String')) returns loc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loc


% --- Executes during object creation, after setting all properties.
function loc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rdrop.
function Rdrop_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选中UMI场景类型的触发事件
if get(handles.Rdrop,'Value')==1 %如果选中random drop
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
else %如果选中manhattan
    if get(handles.UMA2UMI,'Value')
        set(handles.div_angle,'Enable','on');
        set(handles.div_dist,'Enable','on');
        set(handles.div_angle,'String','0');
        set(handles.div_dist,'String','0');
    end
end
% Hints: contents = cellstr(get(hObject,'String')) returns Rdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rdrop


% --- Executes during object creation, after setting all properties.
function Rdrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function femto_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选中InH场景类型的触发事件
% 只返回该值即可
% if get(handles.femto,'Value')==1 %如果选中Home femto
% set(handles.div_angle,'Enable','off');
% set(handles.div_dist,'Enable','off');
% set(handles.div_angle,'String','');
% set(handles.div_dist,'String','');
% else %如果选中manhattan
%     if get(handles.UMA2UMI,'Value')
%         set(handles.div_angle,'Enable','on');
%         set(handles.div_dist,'Enable','on');
%         set(handles.div_angle,'String','0');
%         set(handles.div_dist,'String','0');
%     end
% end
% Hints: contents = cellstr(get(hObject,'String')) returns Rdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rdrop


% --- Executes during object creation, after setting all properties.
function femto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in RMA.
function RMA_Callback(hObject, eventdata, handles)
% hObject    handle to RMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 单击选中InH场景下的触发事件
% 空置其他同构场景的选中状态
set(handles.UMI,'Value',0);
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
%解锁wraparound选项并重置为不采用wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','on');
set(handles.Rdrop,'Enable','off');%锁死UMI场景的类型
set(handles.femto,'Enable','off');%锁死InH场景的类型
set(handles.slot,'String','5');%重置快照数
%重置基站相关参数
set(handles.p_bs1,'String','43');
set(handles.h_bs1,'String','35');
set(handles.bs_max_gain1,'String','8');
%清空并锁死偏置状态的设置
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% % set(handles.loc,'Value',1);
% % set(handles.loc,'Enable','on');
set(handles.res,'String','50'); %重置频率的设置
%重置频率的设置——900里RMA场景有频率要求
set(handles.f1,'String','7');
set(handles.f2,'String','7.2');
%解锁并重置ISD设置
set(handles.isd,'String','1732');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1);%重置handover状态为采用3dB handover
if get(handles.sysnum,'Value')==1 %如果是双系统则解锁共址的选择，设置默认值为共址
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','on');
end
set(handles.shadow_r,'String','0.5');%重置阴影衰落的相关系数为0.5
% 解锁系统类型选择
set(handles.v_type,'Enable','on');
set(handles.a_type,'Enable','on');

% set(handles.sysnum,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of RMA


% --- Executes on button press in UMA.
function UMA_Callback(hObject, eventdata, handles)
% hObject    handle to UMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 单击选中InH场景下的触发事件
% 空置其他同构场景的选中状态
set(handles.UMI,'Value',0);
set(handles.RMA,'Value',0);
set(handles.InH,'Value',0);
%解锁wraparound选项并重置为不采用wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','on');
set(handles.Rdrop,'Enable','off');%锁死UMI场景的类型
set(handles.femto,'Enable','off');%锁死InH场景的类型
set(handles.slot,'String','5');%重置快照数
%重置基站相关参数
set(handles.p_bs1,'String','43');
set(handles.h_bs1,'String','25');
set(handles.bs_max_gain1,'String','8');
%清空并锁死偏置状态的设置——解锁并设置功能转移到共址状态选择的触发事件中
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','5');%重置分辨率
%重置频率的设置
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%解锁并重置ISD设置
set(handles.isd,'String','200');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1); %重置handover状态为不采用3dB handover
if get(handles.sysnum,'Value')==1 %如果是双系统则解锁共址的选择，设置默认值为共址
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','on');
end
set(handles.shadow_r,'String','0.5'); %重置阴影衰落的相关系数为0.5
% 解锁系统类型选择
set(handles.v_type,'Enable','on');
set(handles.a_type,'Enable','on');
% set(handles.sysnum,'Value',1);

% Hint: get(hObject,'Value') returns toggle state of UMA


% --- Executes on selection change in together.
function together_Callback(hObject, eventdata, handles)
% hObject    handle to together (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns together contents as cell array
%        contents{get(hObject,'Value')} returns selected item from together


% --- Executes during object creation, after setting all properties.
function together_CreateFcn(hObject, eventdata, handles)
% hObject    handle to together (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% 选择干扰类型（同步/异步）的触发事件
% 只需传递选中值即可，无其他约束关系
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function div_angle_Callback(hObject, eventdata, handles)
% hObject    handle to div_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统偏置角度的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of div_angle as text
%        str2double(get(hObject,'String')) returns contents of div_angle as a double


% --- Executes during object creation, after setting all properties.
function div_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to div_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function div_dist_Callback(hObject, eventdata, handles)
% hObject    handle to div_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统偏置距离的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of div_dist as text
%        str2double(get(hObject,'String')) returns contents of div_dist as a double


% --- Executes during object creation, after setting all properties.
function div_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to div_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_Callback(hObject, eventdata, handles)
% hObject    handle to res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统分辨率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of res as text
%        str2double(get(hObject,'String')) returns contents of res as a double


% --- Executes during object creation, after setting all properties.
function res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isd_Callback(hObject, eventdata, handles)
% hObject    handle to isd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统站间距的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of isd as text
%        str2double(get(hObject,'String')) returns contents of isd as a double


% --- Executes during object creation, after setting all properties.
function isd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wrap.
function wrap_Callback(hObject, eventdata, handles)
% hObject    handle to wrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统是否采用wraparound的触发事件
% 需要修改分辨率
if get(handles.UMA,'Value') && get(handles.wrap,'Value')==2
    set(handles.res,'String','10');
elseif get(handles.UMA,'Value') && get(handles.wrap,'Value')==1
    set(handles.res,'String','5');
elseif get(handles.RMA,'Value')
    set(handles.res,'String','80');
end
    
    
% Hints: contents = cellstr(get(hObject,'String')) returns wrap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wrap


% --- Executes during object creation, after setting all properties.
function wrap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置victim系统频率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of f1 as text
%        str2double(get(hObject,'String')) returns contents of f1 as a double


% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function band_Callback(hObject, eventdata, handles)
% hObject    handle to band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统带宽的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of band as text
%        str2double(get(hObject,'String')) returns contents of band as a double


% --- Executes during object creation, after setting all properties.
function band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_ue_Callback(hObject, eventdata, handles)
% hObject    handle to p_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统UE发射功率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of p_ue as text
%        str2double(get(hObject,'String')) returns contents of p_ue as a double


% --- Executes during object creation, after setting all properties.
function p_ue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nf_ue_Callback(hObject, eventdata, handles)
% hObject    handle to nf_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统UE侧噪声系数的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of nf_ue as text
%        str2double(get(hObject,'String')) returns contents of nf_ue as a double


% --- Executes during object creation, after setting all properties.
function nf_ue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_type_Callback(hObject, eventdata, handles)
% hObject    handle to v_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择victim系统的类型触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of v_type as text
%        str2double(get(hObject,'String')) returns contents of v_type as a double


% --- Executes during object creation, after setting all properties.
function v_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_bs1_Callback(hObject, eventdata, handles)
% hObject    handle to h_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置victim系统及同构干扰系统的基站高度的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of h_bs1 as text
%        str2double(get(hObject,'String')) returns contents of h_bs1 as a double


% --- Executes during object creation, after setting all properties.
function h_bs1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f2_Callback(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置aggressor系统频率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of f2 as text
%        str2double(get(hObject,'String')) returns contents of f2 as a double


% --- Executes during object creation, after setting all properties.
function f2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_bs1_Callback(hObject, eventdata, handles)
% hObject    handle to p_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置victim系统及同构干扰系统的基站发射功率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of p_bs1 as text
%        str2double(get(hObject,'String')) returns contents of p_bs1 as a double


% --- Executes during object creation, after setting all properties.
function p_bs1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nf_bs_Callback(hObject, eventdata, handles)
% hObject    handle to nf_bs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置基站侧噪声系数的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of nf_bs as text
%        str2double(get(hObject,'String')) returns contents of nf_bs as a double


% --- Executes during object creation, after setting all properties.
function nf_bs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf_bs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_type_Callback(hObject, eventdata, handles)
% hObject    handle to a_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置aggressor系统类型的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of a_type as text
%        str2double(get(hObject,'String')) returns contents of a_type as a double


% --- Executes during object creation, after setting all properties.
function a_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shadow_r_Callback(hObject, eventdata, handles)
% hObject    handle to shadow_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置系统阴影衰落系数的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of shadow_r as text
%        str2double(get(hObject,'String')) returns contents of shadow_r as a double


% --- Executes during object creation, after setting all properties.
function shadow_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shadow_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in handover.
function handover_Callback(hObject, eventdata, handles)
% hObject    handle to handover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择是否采用3dB handover的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: contents = cellstr(get(hObject,'String')) returns handover contents as cell array
%        contents{get(hObject,'Value')} returns selected item from handover


% --- Executes during object creation, after setting all properties.
function handover_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_bs2_Callback(hObject, eventdata, handles)
% hObject    handle to p_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置异构干扰系统的基站发射功率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of p_bs2 as text
%        str2double(get(hObject,'String')) returns contents of p_bs2 as a double


% --- Executes during object creation, after setting all properties.
function p_bs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_bs2_Callback(hObject, eventdata, handles)
% hObject    handle to h_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置aggressor系统基站高度的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of h_bs2 as text
%        str2double(get(hObject,'String')) returns contents of h_bs2 as a double


% --- Executes during object creation, after setting all properties.
function h_bs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in isphoto.
function isphoto_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择是否输出中间结果图像的触发事件
% 只需传递选中值即可，无其他约束关系
% Hint: get(hObject,'Value') returns toggle state of isphoto


% --- Executes on button press in issave.
function issave_Callback(hObject, eventdata, handles)
% hObject    handle to issave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择是否保存波束赋型矩阵的触发事件
% 只需传递选中值即可，无其他约束关系
% Hint: get(hObject,'Value') returns toggle state of issave


% --- Executes on button press in use_cache.
function use_cache_Callback(hObject, eventdata, handles)
% hObject    handle to use_cache (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置是否使用波束赋型缓存的触发事件
% 只需传递选中值即可，无其他约束关系
% Hint: get(hObject,'Value') returns toggle state of use_cache


% --- Executes on button press in ue_bf.
function ue_bf_Callback(hObject, eventdata, handles)
% hObject    handle to ue_bf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择UE是否采用BF触发事件
if get(hObject,'Value') == 1 %不考虑只有UE有波束赋型，而BS无波束赋型的情况。所以UE一旦选择波束赋型，BS也一定选择BF
    set(handles.bs_bf,'Value',1)
end
% Hint: get(hObject,'Value') returns toggle state of ue_bf


% --- Executes on button press in bs_bf.
function bs_bf_Callback(hObject, eventdata, handles)
% hObject    handle to bs_bf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择BS是否采用BF触发事件
if get(hObject,'Value') == 0 %同上，BS一旦没有波束赋型，UE一定没有波束赋型
    set(handles.ue_bf,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of bs_bf



function bs_array_out_Callback(hObject, eventdata, handles)
% hObject    handle to bs_array_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置室外场景BS天线阵的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_array_out as text
%        str2double(get(hObject,'String')) returns contents of bs_array_out as a double


% --- Executes during object creation, after setting all properties.
function bs_array_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_array_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_max_gain_Callback(hObject, eventdata, handles)
% hObject    handle to ue_max_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置UE侧最大单元增益的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of ue_max_gain as text
%        str2double(get(hObject,'String')) returns contents of ue_max_gain as a double


% --- Executes during object creation, after setting all properties.
function ue_max_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_max_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_phi_Callback(hObject, eventdata, handles)
% hObject    handle to bs_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置BS波束水平方向角误差上界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_phi as text
%        str2double(get(hObject,'String')) returns contents of bs_phi as a double


% --- Executes during object creation, after setting all properties.
function bs_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_phi_Callback(hObject, eventdata, handles)
% hObject    handle to ue_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置UE波束水平方向角误差上界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of ue_phi as text
%        str2double(get(hObject,'String')) returns contents of ue_phi as a double


% --- Executes during object creation, after setting all properties.
function ue_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_max_gain2_Callback(hObject, eventdata, handles)
% hObject    handle to bs_max_gain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置异构系统基站最大发射功率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_max_gain2 as text
%        str2double(get(hObject,'String')) returns contents of bs_max_gain2 as a double


% --- Executes during object creation, after setting all properties.
function bs_max_gain2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_max_gain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_max_gain1_Callback(hObject, eventdata, handles)
% hObject    handle to bs_max_gain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置victim系统及同构干扰系统基站发射功率的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_max_gain1 as text
%        str2double(get(hObject,'String')) returns contents of bs_max_gain1 as a double


% --- Executes during object creation, after setting all properties.
function bs_max_gain1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_max_gain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_theta_Callback(hObject, eventdata, handles)
% hObject    handle to bs_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置BS波束垂直角误差上界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_theta as text
%        str2double(get(hObject,'String')) returns contents of bs_theta as a double


% --- Executes during object creation, after setting all properties.
function bs_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_theta_Callback(hObject, eventdata, handles)
% hObject    handle to ue_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置UE波束垂直角误差上界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of ue_theta as text
%        str2double(get(hObject,'String')) returns contents of ue_theta as a double


% --- Executes during object creation, after setting all properties.
function ue_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_lower_Callback(hObject, eventdata, handles)
% hObject    handle to acir_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置研究ACIR下界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of acir_lower as text
%        str2double(get(hObject,'String')) returns contents of acir_lower as a double


% --- Executes during object creation, after setting all properties.
function acir_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_step_Callback(hObject, eventdata, handles)
% hObject    handle to acir_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置ACIR研究不仅的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of acir_step as text
%        str2double(get(hObject,'String')) returns contents of acir_step as a double


% --- Executes during object creation, after setting all properties.
function acir_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slot_Callback(hObject, eventdata, handles)
% hObject    handle to slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置快照数的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of slot as text
%        str2double(get(hObject,'String')) returns contents of slot as a double


% --- Executes during object creation, after setting all properties.
function slot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_upper_Callback(hObject, eventdata, handles)
% hObject    handle to acir_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置所研究的ACIR上界的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of acir_upper as text
%        str2double(get(hObject,'String')) returns contents of acir_upper as a double


% --- Executes during object creation, after setting all properties.
function acir_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_Callback(hObject, eventdata, handles)
% hObject    handle to acir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置ACIR校准值的触发事件
% 只需传递设定值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of acir as text
%        str2double(get(hObject,'String')) returns contents of acir as a double


% --- Executes during object creation, after setting all properties.
function acir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in distribution.
function distribution_Callback(hObject, eventdata, handles)
% hObject    handle to distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择波束误差分布的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: contents = cellstr(get(hObject,'String')) returns distribution contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distribution


% --- Executes during object creation, after setting all properties.
function distribution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in still.
function still_Callback(hObject, eventdata, handles)
% hObject    handle to still (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择UE静止的触发事件

set(handles.move,'Value',0); % 需要空置UE移动的选项
set(handles.speed,'Enable','off'); %锁死速度设置选项
set(handles.TTI,'Enable','off'); %锁死仿真时长设置
% Hint: get(hObject,'Value') returns toggle state of still



function speed_Callback(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置移动速度的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of speed as text
%        str2double(get(hObject,'String')) returns contents of speed as a double


% --- Executes during object creation, after setting all properties.
function speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TTI_Callback(hObject, eventdata, handles)
% hObject    handle to TTI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置仿真时长的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of TTI as text
%        str2double(get(hObject,'String')) returns contents of TTI as a double


% --- Executes during object creation, after setting all properties.
function TTI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TTI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in move.
function move_Callback(hObject, eventdata, handles)
% hObject    handle to move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 选择UE移动的触发事件
set(handles.still,'Value',0); %空置UE静止的选项
set(handles.speed,'Enable','on'); %解锁UE速度选项
set(handles.TTI,'Enable','on'); %解锁仿真时间
% Hint: get(hObject,'Value') returns toggle state of move



function bs_array_in_Callback(hObject, eventdata, handles)
% hObject    handle to bs_array_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置室内基站天线阵规模的触发事件
% 只需传递选中值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of bs_array_in as text
%        str2double(get(hObject,'String')) returns contents of bs_array_in as a double


% --- Executes during object creation, after setting all properties.
function bs_array_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_array_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L_2_F_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 异构场景中后者干扰前者
set(handles.F_2_L,'Value',0); %前者干扰后者与后者干扰前者只能选择其中一个
if get(handles.UMA2UMI,'Value')+get(handles.UMA2InH,'Value')>0 %如果异构系统中有UMA
    set(handles.v_type,'Enable','on');
    set(handles.a_type,'Enable','off');
    set(handles.v_type,'Value',1);
    set(handles.a_type,'Value',1);
end
% Hint: get(hObject,'Value') returns toggle state of L_2_F

function F_2_L_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%异构场景中前者干扰后者
set(handles.L_2_F,'Value',0); %前者干扰后者与后者干扰前者只能选择其中一个
if get(handles.UMA2UMI,'Value')+get(handles.UMA2InH,'Value')>0 %如果异构系统中有UMA
    set(handles.v_type,'Enable','off');
    set(handles.a_type,'Enable','on');
    set(handles.v_type,'Value',1);
    set(handles.a_type,'Value',1);
end
% Hint: get(hObject,'Value') returns toggle state of F_2_L

function UEs_per_slot_Callback(hObject, eventdata, handles)
% hObject    handle to UEs_per_slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 设置每次快照服务的UE数
% 只需传递设定值即可，无其他约束关系
% Hints: get(hObject,'String') returns contents of UEs_per_slot as text
%        str2double(get(hObject,'String')) returns contents of UEs_per_slot as a double


% --- Executes during object creation, after setting all properties.
function UEs_per_slot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UEs_per_slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
