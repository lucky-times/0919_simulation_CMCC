function varargout = LTE_GUI_show_UEs_and_cells(varargin)
% LTE_GUI_SHOW_UES_AND_CELLS MATLAB code for LTE_GUI_show_UEs_and_cells.fig
%      LTE_GUI_SHOW_UES_AND_CELLS, by itself, creates a new LTE_GUI_SHOW_UES_AND_CELLS or raises the existing
%      singleton*.
%
%      H = LTE_GUI_SHOW_UES_AND_CELLS returns the handle to a new LTE_GUI_SHOW_UES_AND_CELLS or the handle to
%      the existing singleton*.
%
%      LTE_GUI_SHOW_UES_AND_CELLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LTE_GUI_SHOW_UES_AND_CELLS.M with the given input arguments.
%
%      LTE_GUI_SHOW_UES_AND_CELLS('Property','Value',...) creates a new LTE_GUI_SHOW_UES_AND_CELLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LTE_GUI_show_UEs_and_cells_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LTE_GUI_show_UEs_and_cells_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LTE_GUI_show_UEs_and_cells

% Last Modified by GUIDE v2.5 17-Jan-2012 19:00:31

% Plots eNodeB and UE positions as read from the laoded simulation results
%
% (c) Josep Colom Ikuno, INTHFT, 2012

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LTE_GUI_show_UEs_and_cells_OpeningFcn, ...
                   'gui_OutputFcn',  @LTE_GUI_show_UEs_and_cells_OutputFcn, ...
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


% --- Executes just before LTE_GUI_show_UEs_and_cells is made visible.
function LTE_GUI_show_UEs_and_cells_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LTE_GUI_show_UEs_and_cells (see VARARGIN)

% Choose default command line output for LTE_GUI_show_UEs_and_cells
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Save user data (the simulation results)
simulation_data = varargin{1};

if length(varargin) > 1
    simulation_data.sim_results_GUI_handle     = varargin{2};
    set(handles.change_results_GUI_checkbox,'Value',true);
else
    simulation_data.sim_results_GUI_handle     = [];
    set(handles.change_results_GUI_checkbox,'Value',false);
end

set(hObject,'UserData',simulation_data);

UE_list = cell(1,length(simulation_data.UEs));
for u_=1:length(simulation_data.UEs)
    UE_list{u_} = sprintf('%g',u_);
end

cell_list = cell(1,length(simulation_data.eNodeB_sites));
for c_=1:length(simulation_data.eNodeB_sectors)
    cell_list{c_} = sprintf('%g',c_);
end

simulation_data.UE_list   = UE_list;
simulation_data.cell_list = cell_list;

% Initialize checkboxes
set(handles.cell_show_cell_areas,'Value',true);
set(handles.cell_all,'Value',true);
set(handles.cell_show_ids,'Value',false);
set(handles.UE_all,'Value',true);

% Fill list boxes
set(handles.UE_listbox,'String',UE_list);
set(handles.UE_listbox,'Max',length(simulation_data.UEs));
set(handles.UE_listbox,'Min',0);

set(handles.cell_listbox,'String',cell_list);
set(handles.cell_listbox,'Max',length(simulation_data.eNodeB_sectors));
set(handles.cell_listbox,'Min',0);

if isfield(simulation_data.SYS_config,'default_shown_GUI_cells') && ~isempty(simulation_data.SYS_config.default_shown_GUI_cells)
    set(handles.cell_listbox,'Value',simulation_data.SYS_config.default_shown_GUI_cells);
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    set(handles.cell_default,'Value',true);
else
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeB_sectors));
    set(handles.UE_listbox,'Value',1:length(simulation_data.UEs));
    set(handles.cell_default,'Value',false);
end

% Plot main plot
plot_cells_and_UEs(handles)

% UIWAIT makes LTE_GUI_show_UEs_and_cells wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = LTE_GUI_show_UEs_and_cells_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in UE_listbox.
function UE_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to UE_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UE_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UE_listbox
set(handles.UE_all,'Value',false);
set(handles.UE_none,'Value',false);
plot_cells_and_UEs(handles);

% --- Executes during object creation, after setting all properties.
function UE_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UE_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cell_listbox.
function cell_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell_listbox
set(handles.cell_all,'Value',false);
set(handles.cell_none,'Value',false);
set(handles.cell_default,'Value',false);
simulation_data = get(handles.figure1,'UserData');
if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces = [simulation_data.simulation_traces.UE_traces];
else
    the_UE_traces = simulation_data.the_UE_traces;
end
UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
set(handles.UE_listbox,'Value',UE_in_this_cells);
plot_cells_and_UEs(handles);

% --- Executes during object creation, after setting all properties.
function cell_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in UE_show_ids.
function UE_show_ids_Callback(hObject, eventdata, handles)
% hObject    handle to UE_show_ids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_show_ids
plot_cells_and_UEs(handles);


% --- Executes on button press in cell_show_ids.
function cell_show_ids_Callback(hObject, eventdata, handles)
% hObject    handle to cell_show_ids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_show_ids
plot_cells_and_UEs(handles);


% --- Executes on button press in UE_all.
function UE_all_Callback(hObject, eventdata, handles)
% hObject    handle to UE_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_all
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.UE_none,'Value',false);
    set(handles.UE_listbox,'Value',1:length(simulation_data.UEs));
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in UE_none.
function UE_none_Callback(hObject, eventdata, handles)
% hObject    handle to UE_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UE_none
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.UE_all,'Value',false);
    set(handles.UE_listbox,'Value',[]);
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in cell_all.
function cell_all_Callback(hObject, eventdata, handles)
% hObject    handle to cell_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_all
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_none,'Value',false);
    set(handles.cell_default,'Value',false);
    set(handles.cell_listbox,'Value',1:length(simulation_data.eNodeB_sectors));
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    plot_cells_and_UEs(handles);
end


% --- Executes on button press in cell_none.
function cell_none_Callback(hObject, eventdata, handles)
% hObject    handle to cell_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_none
simulation_data = get(handles.figure1,'UserData');
if get(hObject,'Value')
    set(handles.cell_all,'Value',false);
    set(handles.cell_default,'Value',false);
    set(handles.cell_listbox,'Value',[]);
    if isfield(simulation_data,'simulation_traces')
        % non-compact results
        the_UE_traces = [simulation_data.simulation_traces.UE_traces];
    else
        the_UE_traces = simulation_data.the_UE_traces;
    end
    UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
    set(handles.UE_listbox,'Value',UE_in_this_cells);
    plot_cells_and_UEs(handles);
end

function plot_cells_and_UEs(handles,varargin)

if isempty(varargin)
    plot_in_GUI = true;
else
    plot_in_GUI = varargin{1};
end

% Load data
simulation_data    = get(handles.figure1,'UserData');
eNodeB_sites            = simulation_data.eNodeB_sites;
eNodeB_sectors    = simulation_data.eNodeB_sectors;
UEs                = simulation_data.UEs;
networkPathlossMap = simulation_data.networkPathlossMap;

UEs_to_plot_lin    = get(handles.UE_listbox,'Value');
cells_to_plot      = get(handles.cell_listbox,'Value');
plot_cell_ids      = get(handles.cell_show_ids,'Value');
plot_cell_area     = get(handles.cell_show_cell_areas,'Value');
plot_UE_ids        = get(handles.UE_show_ids,'Value');
change_results_GUI = get(handles.change_results_GUI_checkbox,'Value');

if isfield(simulation_data,'simulation_traces')
    % non-compact results
    the_UE_traces = [simulation_data.simulation_traces.UE_traces];
else
    the_UE_traces = simulation_data.the_UE_traces;
end

if change_results_GUI && ~isempty(simulation_data.sim_results_GUI_handle) && ishandle(simulation_data.sim_results_GUI_handle)
    sim_results_GUI_handle_guidata = guidata(simulation_data.sim_results_GUI_handle);
    set(sim_results_GUI_handle_guidata.cell_listbox,'Value',cells_to_plot);
    set(sim_results_GUI_handle_guidata.cell_all,'Value',false);
    set(sim_results_GUI_handle_guidata.cell_none,'Value',false);
end

UE_pos             = reshape([UEs.pos],2,length(UEs))';
site_pos           = reshape([eNodeB_sites.pos],2,length(eNodeB_sites))';

if plot_in_GUI
    the_axes = handles.map_axes;
else
    the_new_figure = figure;
    the_axes       = axes('Parent',the_new_figure);
end

cla(the_axes);
hold(the_axes,'all');
set(the_axes,'YDir','normal');

% Shade the cell area of the specified cell/s
shaded_DotColor = 0.85*[1 1 1];
if plot_cell_area
    pos_to_shade = false(size(networkPathlossMap.sector_assignment));
    for c_ = cells_to_plot
        if simulation_data.SYS_config.isDouble
            n_sector = networkPathlossMap.num_first_sectors;
        else
            n_sector = length(eNodeB_sectors);
        end
        if c_<=n_sector
            pos_to_shade = pos_to_shade | (networkPathlossMap.sector_assignment==c_);
        else
            pos_to_shade = pos_to_shade | (networkPathlossMap.sector_assignment_double==c_);
        end
    end
    [row,col] = find(pos_to_shade);
    shaded_points = NR_common_pixel_to_pos([col row],networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
    scatter(the_axes,shaded_points(:,1),shaded_points(:,2),'Marker','.','MarkerFaceColor',shaded_DotColor,'MarkerEdgeColor',shaded_DotColor);
end

% Add the position of each UE

UE_DotColor_blue = [0 0 1];
dark_blue        = [0.7 0.7 1];
UE_DotColor_red  = [1 0 0];
dark_red         = [1 0.7 0.7];
UE_DotColor_yellow = [0 1 0];
dark_yellow = [0.7 1 0.7];

if simulation_data.SYS_config.isDouble
    num_first_UEs = networkPathlossMap.num_first_UEs;
    index_UEs1 = 1:num_first_UEs;
    index_UEs2 = num_first_UEs+1:length(UEs);
    UEs_to_plot1 = false(length(UEs),1);
    UEs_to_plot2 = false(length(UEs),1);
    rest_UEs1 = true(length(UEs),1);
    rest_UEs2 = true(length(UEs),1);
    UEs_to_plot_lin1 = intersect(UEs_to_plot_lin,index_UEs1);
    UEs_to_plot_lin2 = intersect(UEs_to_plot_lin,index_UEs2);
    UEs_to_plot1(UEs_to_plot_lin1) = true;
    UEs_to_plot2(UEs_to_plot_lin2) = true;
    rest_UEs1(UEs_to_plot1) = false;
    rest_UEs2(UEs_to_plot2) = false;
    rest_UEs1(index_UEs2) = false;
    rest_UEs2(index_UEs1) = false;
else
    UEs_to_plot = false(length(UEs),1);
    rest_UEs = true(length(UEs),1);
    UEs_to_plot(UEs_to_plot_lin) = true;
    rest_UEs(UEs_to_plot) = false;
end

if simulation_data.SYS_config.isDouble
    scatter(the_axes,UE_pos(rest_UEs1,1),UE_pos(rest_UEs1,2),'Marker','.','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_blue);
    scatter(the_axes,UE_pos(rest_UEs2,1),UE_pos(rest_UEs2,2),'Marker','.','MarkerFaceColor',dark_yellow,'MarkerEdgeColor',dark_yellow);
    scatter(the_axes,UE_pos(UEs_to_plot1,1),UE_pos(UEs_to_plot1,2),'Marker','.','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_blue);
    scatter(the_axes,UE_pos(UEs_to_plot2,1),UE_pos(UEs_to_plot2,2),'Marker','.','MarkerFaceColor',UE_DotColor_yellow,'MarkerEdgeColor',UE_DotColor_yellow);
else
    scatter(the_axes,UE_pos(rest_UEs,1),UE_pos(rest_UEs,2),'Marker','.','MarkerFaceColor',dark_blue,'MarkerEdgeColor',dark_blue);
    scatter(the_axes,UE_pos(UEs_to_plot,1),UE_pos(UEs_to_plot,2),'Marker','.','MarkerFaceColor',UE_DotColor_blue,'MarkerEdgeColor',UE_DotColor_blue);
end

if plot_UE_ids
    for u_=UEs_to_plot_lin
        text(UEs(u_).pos(1)+15*1,UEs(u_).pos(2),num2str(UEs(u_).id),'FontSize',8);
    end
end

% Plot a line that tells where the antennas are pointing
antenna_LineColor = 'blue';
for b_=1:length(eNodeB_sites)
    origin = eNodeB_sites(b_).pos;
    for s_=1:length(eNodeB_sites(b_).sectors)
        switch eNodeB_sites(b_).site_type
            case 'macro'
                angle = wrapTo360(-eNodeB_sites(b_).sectors(s_).azimuth+90);
                vector_length = simulation_data.SYS_config.ISD * 0.1;
            case 'micro'
                angle = eNodeB_sites(b_).sectors(s_).azimuth;
                if ~simulation_data.SYS_config.isManhattan
                    vector_length = simulation_data.SYS_config.ISD * 0.1;
                else
                    vector_length = 0;
                end
            case 'indoor'
                angle = wrapTo360(-eNodeB_sites(b_).sectors(s_).azimuth+90);
                vector_length = 0;
        end
        vector = vector_length*[ cosd(angle) sind(angle) ];
        destiny = vector + origin;
        plot([origin(1) destiny(1)],[origin(2) destiny(2)],antenna_LineColor);
    end
end

% Plot all of the sites' positions
site_MarkerEdgeColor = 'black';
site_MarkerFaceColor = 'red';
site_MarkerFaceColor2 = 'yellow';

if simulation_data.SYS_config.isDouble
    n_sites = networkPathlossMap.num_first_sites;
else
    n_sites = length(eNodeB_sites);
end
for b_ = 1:n_sites
    scatter(eNodeB_sites(b_).pos(1),eNodeB_sites(b_).pos(2),'MarkerEdgeColor',site_MarkerEdgeColor,'MarkerFaceColor',site_MarkerFaceColor);
    text(eNodeB_sites(b_).pos(1),eNodeB_sites(b_).pos(2),num2str(b_),'Color','k');
end
if simulation_data.SYS_config.isDouble
    for b_ = n_sites+1:length(eNodeB_sites)
        scatter(eNodeB_sites(b_).pos(1),eNodeB_sites(b_).pos(2),'MarkerEdgeColor',site_MarkerEdgeColor,'MarkerFaceColor',site_MarkerFaceColor2);
        text(eNodeB_sites(b_).pos(1),eNodeB_sites(b_).pos(2),num2str(b_),'Color','k');
    end
end

% Mark with a text the center of the cells
if plot_cell_ids
    if simulation_data.SYS_config.isDouble
        n_sector = networkPathlossMap.num_first_sectors;
    else
        n_sector = length(eNodeB_sectors);
    end
    for c_ = 1:n_sector % cells_to_plot
        text(networkPathlossMap.sector_centers(c_,1),networkPathlossMap.sector_centers(c_,2),num2str(c_),'HorizontalAlignment','center','Verticalalignment','middle','Color',0.50*[1 1 1]);
    end
    if simulation_data.SYS_config.isDouble
        for c_ = n_sector+1:length(eNodeB_sectors)
            text(networkPathlossMap.sector_centers_double(c_-n_sector,1),networkPathlossMap.sector_centers_double(c_-n_sector,2),num2str(c_),'HorizontalAlignment','center','Verticalalignment','middle','Color',0.50*[1 1 1]);
        end
    end
end

% 画网格
hex_color = 'red';
hex_color2 = 'yellow';
rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)];

if simulation_data.SYS_config.isDouble
    n_sites = networkPathlossMap.num_first_sites;
else
    n_sites = length(eNodeB_sites);
end
switch simulation_data.SYS_config.scene_type
    case {'UMA','RMa','UMI','UMa_to_UMi'}
        ISD = simulation_data.SYS_config.ISD;
        if ~simulation_data.SYS_config.isManhattan
            for b_ = 1:n_sites
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    for i_ = 1:7
                        tmp_hex(i_,:) = eNodeB_sites(b_).sector_centre(s_).pos + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % 乘[1;0]表示只取第一列
                    end
                    plot(tmp_hex(:,1),tmp_hex(:,2),hex_color);
                end
            end
            if simulation_data.SYS_config.isDouble
                for b_ = n_sites+1:length(eNodeB_sites)
                    for s_ = 1:length(eNodeB_sites(b_).sectors)
                        for i_ = 1:7
                            tmp_hex(i_,:) = eNodeB_sites(b_).sector_centre(s_).pos + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % 乘[1;0]表示只取第一列
                        end
                        plot(tmp_hex(:,1),tmp_hex(:,2),hex_color2);
                    end
                end
            end
        else
            switch simulation_data.SYS_config.scene_type
                case 'UMa_to_UMi'
                    for b_ = n_sites+1:length(eNodeB_sites)
                        for s_ = 1:length(eNodeB_sites(b_).sectors)
                            for i_ = 1:7
                                tmp_hex(i_,:) = eNodeB_sites(b_).sector_centre(s_).pos + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % 乘[1;0]表示只取第一列
                            end
                            plot(tmp_hex(:,1),tmp_hex(:,2),hex_color2);
                        end
                    end
            end
        end
    case {'UMa_to_InH','UMi_to_InH'}
        ISD = simulation_data.SYS_config.ISD;
        if simulation_data.SYS_config.isFemto
            for b_ = 1:n_sites
                femto_room = eNodeB_sites(b_).femto_room;
                room_x = [femto_room(1,1),femto_room(1,1),femto_room(1,2),femto_room(1,2),femto_room(1,1)];
                room_y = [femto_room(2,1),femto_room(2,2),femto_room(2,2),femto_room(2,1),femto_room(2,1)];
                plot(room_x,room_y,hex_color);
            end
        end
        for b_ = n_sites+1:length(eNodeB_sites)
            for s_ = 1:length(eNodeB_sites(b_).sectors)
                for i_ = 1:7
                    tmp_hex(i_,:) = eNodeB_sites(b_).sector_centre(s_).pos + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % 乘[1;0]表示只取第一列
                end
                plot(tmp_hex(:,1),tmp_hex(:,2),hex_color2);
            end
        end
    case 'InH'
        if simulation_data.SYS_config.isFemto
            for b_ = 1:n_sites
                femto_room = eNodeB_sites(b_).femto_room;
                room_x = [femto_room(1,1),femto_room(1,1),femto_room(1,2),femto_room(1,2),femto_room(1,1)];
                room_y = [femto_room(2,1),femto_room(2,2),femto_room(2,2),femto_room(2,1),femto_room(2,1)];
                plot(room_x,room_y,hex_color);
            end
            if simulation_data.SYS_config.isDouble
                for b_ = n_sites+1:length(eNodeB_sites)
                    femto_room = eNodeB_sites(b_).femto_room;
                    room_x = [femto_room(1,1),femto_room(1,1),femto_room(1,2),femto_room(1,2),femto_room(1,1)];
                    room_y = [femto_room(2,1),femto_room(2,2),femto_room(2,2),femto_room(2,1),femto_room(2,1)];
                    plot(room_x,room_y,hex_color2);
                end
            end
        end
end

xlabel(the_axes,'x pos [m]');
ylabel(the_axes,'y pos [m]');
xlim(the_axes,networkPathlossMap.roi_x);
ylim(the_axes,networkPathlossMap.roi_y);
title(the_axes,sprintf('eNodeB and UE positions'));

if ~plot_in_GUI
    axis(the_axes,'equal');
end

% Release the hold on the axes
hold(the_axes,'off');



% --- Executes on button press in cell_show_cell_areas.
function cell_show_cell_areas_Callback(hObject, eventdata, handles)
% hObject    handle to cell_show_cell_areas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_show_cell_areas
plot_cells_and_UEs(handles);


% --- Executes on button press in open_plot_in_new_figure.
function open_plot_in_new_figure_Callback(hObject, eventdata, handles)
% hObject    handle to open_plot_in_new_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_cells_and_UEs(handles,false);


% --- Executes on button press in change_results_GUI_checkbox.
function change_results_GUI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to change_results_GUI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of change_results_GUI_checkbox
simulation_data    = get(handles.figure1,'UserData');
current_status     = get(hObject,'Value');
if ~isempty(simulation_data.sim_results_GUI_handle) && ishandle(simulation_data.sim_results_GUI_handle)
    if current_status
        cells_to_plot = get(handles.cell_listbox,'Value');
        sim_results_GUI_handle_guidata = guidata(simulation_data.sim_results_GUI_handle);
        set(sim_results_GUI_handle_guidata.cell_listbox,'Value',cells_to_plot);
        set(sim_results_GUI_handle_guidata.cell_all,'Value',false);
        set(sim_results_GUI_handle_guidata.cell_none,'Value',false);
    end
else
    set(handles.change_results_GUI_checkbox,'Value',~current_status);
end


% --- Executes on button press in cell_default.
function cell_default_Callback(hObject, eventdata, handles)
% hObject    handle to cell_default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cell_default
simulation_data = get(handles.figure1,'UserData');

if isfield(simulation_data.SYS_config,'default_shown_GUI_cells') && ~isempty(simulation_data.SYS_config.default_shown_GUI_cells)
    if get(hObject,'Value')
        set(handles.cell_listbox,'Value',simulation_data.SYS_config.default_shown_GUI_cells);
        if isfield(simulation_data,'simulation_traces')
            % non-compact results
            the_UE_traces = [simulation_data.simulation_traces.UE_traces];
        else
            the_UE_traces = simulation_data.the_UE_traces;
        end
        UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
        set(handles.UE_listbox,'Value',UE_in_this_cells);
        set(handles.cell_all,'Value',false);
        set(handles.cell_none,'Value',false);
        
        plot_cells_and_UEs(handles);
    end
else
    set(handles.cell_default,'Value',~get(hObject,'Value'));
end


% if get(hObject,'Value')
%     set(handles.cell_all,'Value',false);
%     set(handles.cell_default,'Value',false);
%     set(handles.cell_listbox,'Value',[]);
%     if isfield(simulation_data,'simulation_traces')
%         % non-compact results
%         the_UE_traces = [simulation_data.simulation_traces.UE_traces];
%     else
%         the_UE_traces = simulation_data.the_UE_traces;
%     end
%     UE_in_this_cells = find(utils.miscUtils.get_UEs_in_given_cells(get(handles.cell_listbox,'Value'),the_UE_traces));
%     set(handles.UE_listbox,'Value',UE_in_this_cells);
%     plot_cells_and_UEs(handles);
% end
