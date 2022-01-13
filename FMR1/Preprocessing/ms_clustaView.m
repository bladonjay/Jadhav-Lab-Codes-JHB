function varargout = ms_clustaView(varargin)
% MS_CLUSTAVIEW MATLAB code for ms_clustaView.fig
%      MS_CLUSTAVIEW, by itself, creates a new MS_CLUSTAVIEW or raises the existing
%      singleton*.
%
%      H = MS_CLUSTAVIEW returns the handle to a new MS_CLUSTAVIEW or the handle to
%      the existing singleton*.
%
%      MS_CLUSTAVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MS_CLUSTAVIEW.M with the given input arguments.
%
%      MS_CLUSTAVIEW('Property','Value',...) creates a new MS_CLUSTAVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ms_clustaView_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ms_clustaView_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_clustaView

% Last Modified by GUIDE v2.5 15-Nov-2018 15:48:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ms_clustaView_OpeningFcn, ...
                   'gui_OutputFcn',  @ms_clustaView_OutputFcn, ...
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


% --- Executes just before ms_clustaView is made visible.
function ms_clustaView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_clustaView (see VARARGIN)

% Choose default command line output for ms_clustaView
handles.output = hObject;
setappdata(handles.figure1,'clustData',[])
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_clustaView wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ms_clustaView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in x_pop.
function x_pop_Callback(hObject, eventdata, handles)
% hObject    handle to x_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_plot(handles)


% --- Executes during object creation, after setting all properties.
function x_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in y_pop.
function y_pop_Callback(hObject, eventdata, handles)
% hObject    handle to y_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_plot(handles)


% --- Executes during object creation, after setting all properties.
function y_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function open_menu_Callback(hObject, eventdata, handles)
% hObject    handle to open_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn,fileDir] = uigetfile('*.mat','Choose mountainsort-derived matclust params file');
if isempty(fn)
    return;
end
clustDat = load([fileDir filesep fn]);
if ~isfield(clustDat,'filedata')
    errordlg('Invalid file chosen')
    return;
end
fd = clustDat.filedata;
popStr = fd.paramnames;
set(handles.y_pop,'String',popStr,'Value',1)
set(handles.x_pop,'String',popStr,'Value',2)
cidx = find(strcmpi(fd.paramnames,'Cluster'));
if ~isempty(cidx)
    clusts = sort(unique(fd.params(:,cidx)));
    tblDat = num2cell([clusts logical(zeros(size(clusts,1),1))]);
else
    tblDat = {1,logical(0)};
end
for k=1:size(tblDat,1)
    tblDat{k,2} = logical(tblDat{k,2});
end
set(handles.clust_list,'Data',tblDat)
colors = lines(size(tblDat,1));
setappdata(handles.figure1,'plotColors',colors)
setappdata(handles.figure1,'clustData',fd);
update_plot(handles)



% --- Executes when entered data in editable cell(s) in clust_list.
function clust_list_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to clust_list (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
update_plot(handles)


function update_plot(handles)
clustDat = getappdata(handles.figure1,'clustData');
if isempty(clustDat)
    return;
end
yPopStr = get(handles.y_pop,'String');
yVal = get(handles.y_pop,'Value');
y_pop = yPopStr{yVal};
xPopStr = get(handles.x_pop,'String');
xVal = get(handles.x_pop,'Value');
x_pop = xPopStr{xVal};

paramnames = clustDat.paramnames;
yidx = strcmpi(paramnames,y_pop);
xidx = strcmpi(paramnames,x_pop);

ydat = clustDat.params(:,yidx);
xdat = clustDat.params(:,xidx);

cidx = find(strcmpi(paramnames,'Cluster'));
if ~isempty(cidx)
    cdat = clustDat.params(:,cidx);
else
    cdat = ones(size(ydat,1),1);
end
tblDat = get(handles.clust_list,'Data');
colors = getappdata(handles.figure1,'plotColors');
cla(handles.plot_ax)
hold(handles.plot_ax,'on')
for k=1:size(tblDat,1)
    clust = tblDat{k,1};
    if tblDat{k,2}
        continue;
    end
    idx = cdat==clust;
    
    plot(xdat(idx),ydat(idx),'.','Color',colors(k,:))
end
