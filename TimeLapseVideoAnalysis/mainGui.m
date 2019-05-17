function varargout = mainGui(varargin)
    % MAINGUI MATLAB code for mainGui.fig
    %      MAINGUI, by itself, creates a new MAINGUI or raises the existing
    %      singleton*.
    %
    %      H = MAINGUI returns the handle to a new MAINGUI or the handle to
    %      the existing singleton*.
    %
    %      MAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in MAINGUI.M with the given input arguments.
    %
    %      MAINGUI('Property','Value',...) creates a new MAINGUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before mainGui_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to mainGui_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help mainGui

    % Last Modified by GUIDE v2.5 12-May-2017 23:24:21

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @mainGui_OpeningFcn, ...
                       'gui_OutputFcn',  @mainGui_OutputFcn, ...
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
end

% --- Executes just before mainGui is made visible.
function mainGui_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to mainGui (see VARARGIN)

    % Choose default command line output for mainGui
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes mainGui wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = mainGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end

function step3_executeButton_Callback(hObject, eventdata, handles)
    GUI_2_handle = displayResults_step3; % open GUI_2 and save the handle
    %CurButtonGroup = get(hObject, 'parent');
    %mainWindow = get(CurButtonGroup, 'parent');
    %delete(get(mainWindow, 'parent')); % close GUI_1 
end

function step2_executeButton_Callback(hObject, eventdata, handles)
    stackToCsv_step2();
end

function step1_executeButton_Callback(hObject, eventdata, handles)
    analyseAndTrackCells_step1()
end
