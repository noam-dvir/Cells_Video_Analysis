
function varargout = displayResults_step3(varargin)
    % DISPLAYRESULTS_STEP3 MATLAB code for displayResults_step3.fig
    %      DISPLAYRESULTS_STEP3, by itself, creates a new DISPLAYRESULTS_STEP3 or raises the existing
    %      singleton*.
    %
    %      H = DISPLAYRESULTS_STEP3 returns the handle to a new DISPLAYRESULTS_STEP3 or the handle to
    %      the existing singleton*.
    %
    %      DISPLAYRESULTS_STEP3('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in DISPLAYRESULTS_STEP3.M with the given input arguments.
    %
    %      DISPLAYRESULTS_STEP3('Property','Value',...) creates a new DISPLAYRESULTS_STEP3 or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before displayResults_step3_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to displayResults_step3_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help displayResults_step3

    % Last Modified by GUIDE v2.5 15-May-2017 15:41:10

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @displayResults_step3_OpeningFcn, ...
                       'gui_OutputFcn',  @displayResults_step3_OutputFcn, ...
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

% --- Executes just before displayResults_step3 is made visible.
function displayResults_step3_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to displayResults_step3 (see VARARGIN)

    % Choose default command line output for displayResults_step3
    handles.output = hObject;
    handles.currRatio = 0;

    while 1
        [movieFile, movieBaseFolder] = uigetfile('*.tif','select the segmented movie to analyes');
        
        % only for debug
        %movieBaseFolder = '/Users/eyal/Documents/Dropbox/CS/0_My_Semester/FinalProj/myProj/test/GUI2/data/dataset7/';
        %movieFile = 'out-segmentedMovie.tif';
        if exist(strcat(movieBaseFolder, movieFile), 'file')
            break;
        end
        uiwait(msgbox('input movie was not selected', 'Error','error'));
    end
    
    handles.baseFolder = movieBaseFolder;
    handles.trackFile = movieFile;
    handles.fullMovieFile = strcat(handles.baseFolder, handles.trackFile);
    
    while 1
        [dbFile, dbBasePFolder] = uigetfile('*.csv','select the .csv file');
        
        % only for debug
        %dbBasePFolder = '/Users/eyal/Documents/Dropbox/CS/0_My_Semester/FinalProj/myProj/test/GUI2/data/dataset7/';
        %dbFile = 'csvCell.csv';
        if exist(strcat(dbBasePFolder, dbFile), 'file')
            break;
        end
        uiwait(msgbox('input DB was not selected', 'Error','error'));
    end
    
    handles.cellMatrix = parseCsv2CellMatrix(strcat(dbBasePFolder,dbFile));
    handles.sharedData = calcSharedData(handles.cellMatrix);
    handles.errorMsg = 'aa';
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes displayResults_step3 wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = displayResults_step3_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

function text_input_Callback(hObject, eventdata, handles)
end

function text_input_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end








%% start from here

function ratio_cellArea_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_cellPerimeter_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_cellPosition_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_cellNeighboursCounter_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_3plot_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_sharedData_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles)
end

function ratio_edgeSize_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_4plot_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_colorSingle_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_color_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_t1_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_menualCorrection_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_chooseLbl_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function ratio_fixDrift_Callback(hObject, eventdata, handles)
    handle_all_ratio(hObject, eventdata, handles);
end

function handle_all_ratio(hObject, eventdata, handles)
    if (handles.currRatio ~= 0)
        set(handles.currRatio, 'Value', 0);
    end
    handles.currRatio = hObject;
    guidata(hObject,handles)
end

function pushbutton1_Callback(hObject, eventdata, handles)
    result = validate_input(handles);
    if result.isValid
        result = callFuncs(hObject, handles);
        uiwait(msgbox(result.msg, 'Cmd Summery'));
    else
        uiwait(msgbox(result.msg, 'Cmd Summery', 'Error', 'error'));
    end
end

function y = validate_input(handles)
    % each cmd has diffrect restrictions for input- this function validate
    % the input for each cmd
    lblVec = str2num(get(handles.text_input, 'String'));
    y.isValid = 1;
    y.msg = 'input valid';
    
    curState = get(handles.currRatio, 'tag');
    maxLbl = size(handles.cellMatrix, 1);
    switch curState
        case 'ratio_sharedData'
        case 'ratio_3plot'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_cellArea'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_cellPerimeter'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_cellPosition'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_cellNeighboursCounter'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_4plot'
            y = validateLblVecExact(lblVec, 1, maxLbl);
        case 'ratio_edgeSize'
            y = validateLblVecExact(lblVec, 2, maxLbl);
        case 'ratio_t1'
        case 'ratio_colorSingle'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_color'
        case 'ratio_menualCorrection'
            y = validateLblVecAbove(lblVec, 0, maxLbl);
        case 'ratio_chooseLbl'
            y = validateFrameInRange(handles, lblVec);
        case 'ratio_fixDrift'
        otherwise
    end
end
% current cmd is stored in handles.currRatio, please see handle_all_ratio()
function y = callFuncs(hObject, handles)
    lblVec = str2num(get(handles.text_input, 'String'));
    
    curState = get(handles.currRatio, 'tag');
    switch curState
        case 'ratio_sharedData'
            plotSharedData(handles.sharedData);
        case 'ratio_3plot'
            plotCellAreaPerimeterNeighbor(handles.cellMatrix, lblVec);
        case 'ratio_cellArea'
            plotCellArea(handles.cellMatrix, lblVec);
        case 'ratio_cellPerimeter'
            plotCellPerimeter(handles.cellMatrix, lblVec);
        case 'ratio_cellPosition'
            plotCellPosition(handles.cellMatrix, handles.sharedData, lblVec);
        case 'ratio_cellNeighboursCounter'
            plotCellNeighbourCounter(handles.cellMatrix, lblVec);
        case 'ratio_4plot'
            plotSingleCellAreaPerimeterPositionNeighbor(handles.cellMatrix, handles.sharedData, lblVec);
        case 'ratio_edgeSize'
            plotAdjacentCellsEdgeSize(handles.cellMatrix, lblVec);
        case 'ratio_t1'
            lookForT1Transitions(handles.cellMatrix, handles);
        case 'ratio_colorSingle'
            plotCellOnStack(lblVec, handles);
        case 'ratio_color'
            colorStackRandom(handles);
        case 'ratio_menualCorrection'
            menualCorrection(handles, lblVec);
        case 'ratio_chooseLbl'
            chooseLbl(hObject, handles, lblVec);
        case 'ratio_fixDrift'
            fixDrift(handles);
        otherwise
    end
    
    y.isValid = 1;
    y.msg = strcat(curState,' was finished');
end

function y = validateLblVecExact(lblVec, i, maxLbl)
    if size(lblVec, 1) == 0 || size(lblVec, 2) == 0
        y.isValid = 0;
        y.msg = strcat('input is not valide- input should contains only numbers');
    elseif size(lblVec, 2) ~= i
        y.isValid = 0;
        y.msg = strcat('input is not valide- number of lbls should be exectly ', num2str(i));
    elseif maxLbl < max(lblVec)
        y.isValid = 0;
        y.msg = 'input is not valide- lbl exceed max lbl from segmentation, peak lower lbl';
    else
        y.isValid = 1;
        y.msg = 'inpus is valid';
    end
end

function y = validateFrameInRange(handles, frame)
    stackInfo = imfinfo(handles.fullMovieFile);
    nFrame = numel(stackInfo);
    
    if size(frame, 1) == 0 || size(frame, 2) == 0
        y.isValid = 0;
        y.msg = strcat('input is not valide- input should contains only numbers');
    elseif size(frame, 2) ~= 1
        y.isValid = 0;
        y.msg = strcat('input is not valide- number of frames should be exectly 1');
    elseif nFrame < frame
        y.isValid = 0;
        y.msg = 'input is not valide- selected frame is out of range, peak lower frame';
    else
        y.isValid = 1;
        y.msg = 'inpus is valid';
    end
end

function y = validateLblVecAbove(lblVec, i, maxLbl)
    y.isValid = 1;
    y.msg = 'eyal';
    if size(lblVec, 2) <= i
        y.isValid = 0;
        y.msg = strcat('input is not valide- number of lbls should be more than- ', num2str(i));
    elseif maxLbl < max(lblVec)
        y.isValid = 0;
        y.msg = 'input is not valide- lbl exceed max lbl from segmentation, peak lower lbl';
    end
end














% this function take the DB and convert it to matrix of MyCell
function y = parseCsv2CellMatrix(csvFile)
    cellMatrix = {};
    M = csvread(csvFile);
    for i = 1: size(M, 1)
        c = MyCell;
        cellData = M(i,:);
        fromMatrix(c, cellData);
        cellMatrix{c.lbl, c.frame} = c;
    end
    y = cellMatrix;
end


function y = calcSharedData(cellMatrix, sharedData)
    nFrame = size(cellMatrix, 2);
    
    cellsPerFrameVec = zeros(1, nFrame);
    
    areaPerFrameVec = zeros(1, nFrame);
    perimeterPerFrameVec = zeros(1, nFrame);
    neighbourPerFrameVec = zeros(1, nFrame);
    
    driftByFrameDataVer1 = zeros(2, nFrame);
    driftByFrameDataVer2 = zeros(2, nFrame); 
    
    for iFrame = 1:nFrame
        for lbl = 1:size(cellMatrix, 1)
            curCell = cellMatrix{lbl, iFrame};
            if (isempty(curCell))
                continue;
            end
            cellsPerFrameVec(iFrame) = cellsPerFrameVec(iFrame) + 1;
            
            areaPerFrameVec(iFrame) = areaPerFrameVec(iFrame) + curCell.area;
            perimeterPerFrameVec(iFrame) = perimeterPerFrameVec(iFrame) + curCell.perimeter;
            neighbourPerFrameVec(iFrame) = neighbourPerFrameVec(iFrame) + getNeighborsSum(curCell);

            driftByFrameDataVer1(:, iFrame) = driftByFrameDataVer1(:, iFrame) + getDriftCord(curCell).';
            driftByFrameDataVer2(:, iFrame) = driftByFrameDataVer2(:, iFrame) + curCell.centroid.';
        end
    end
    
    sharedData.cellsPerFrame = cellsPerFrameVec;
    
    sharedData.areaPerFrame = areaPerFrameVec;
    sharedData.perimeterPerFrame = perimeterPerFrameVec;
    sharedData.neighbourPerFrame = neighbourPerFrameVec;
    % plesae note that for drift we will use sharedData.driftPerFrameVer1
    sharedData.driftPerFrameVer1 = driftByFrameDataVer1(:, :)./[cellsPerFrameVec ; cellsPerFrameVec];
    
    % second drift index- this is not too accurate
    % todo- can be deleted or leave here for future improvments
    meanDriftByFrameDataVer2(:, :) = driftByFrameDataVer2(:, :)./cellsPerFrameVec;
    secondDriftIndex = zeros(2, nFrame);
    for iFrame = 2:nFrame
        secondDriftIndex(:, iFrame) = (meanDriftByFrameDataVer2(:, iFrame) - meanDriftByFrameDataVer2(:, iFrame-1));
    end
    sharedData.driftPerFrameVer2 = secondDriftIndex;
    
    y = sharedData;
end

function y = getAvarageCellArea(sharedData)
    y = sum(sharedData.areaPerFrame)/sum(sharedData.cellsPerFrame);
end

function y = getAvarageCellAreaPerFrame(sharedData)
    y = sharedData.areaPerFrame./sharedData.cellsPerFrame;
end

function y = getAvarageCellPerimeter(sharedData)
    y = sum(sharedData.perimeterPerFrame)/sum(sharedData.cellsPerFrame);
end

function y = getAvarageCellPerimeterPerFrame(sharedData)
    y = sharedData.perimeterPerFrame./sharedData.cellsPerFrame;
end

function y = getAvarageNeighbourCounter(sharedData)
    y = sum(sharedData.neighbourPerFrame)/sum(sharedData.cellsPerFrame);
end

function y = getAvarageNeighbourCounterPerFrame(sharedData)
    y = sharedData.neighbourPerFrame./sharedData.cellsPerFrame;
end


function y = plotSharedData(sharedData)
    sharedDataFigure = figure;  % this figure is the current active automatically
    
    avarageAreaVec = getAvarageCellAreaPerFrame(sharedData);
    avaragePerimeter = getAvarageCellPerimeterPerFrame(sharedData);
    avarageNeighbour =  getAvarageNeighbourCounterPerFrame(sharedData);
    
    nFrames = size(avarageAreaVec, 2);
    subplot(3,1,1);
    plot(1:nFrames, avarageAreaVec);
    title('Avarage Cell Area by Frame');
    
    subplot(3,1,2);
    plot(1:nFrames, avaragePerimeter);
    title('Avarage Cell Perimeter by Frame');
    
    subplot(3,1,3);
    plot(1:nFrames, avarageNeighbour);
    ylim([0 10]);
    title('Avarage # of Neighbour by Frame');
end

function y = plotCellArea(cellMatrix, lblVec)
    areaFigure = figure;  % this figure is the current active automatically
    areaMatrix = getCellAreaMatrix(cellMatrix, lblVec);
    nFrames = size(areaMatrix, 2);
    plot(1:nFrames, areaMatrix);
    title('Area by Frame');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
end

function y = getCellAreaMatrix(cellMatrix, lblVec)
    assert(max(lblVec) <= size(cellMatrix, 1));     % make sure all the lbls from the usee are valid
    nFrames = size(cellMatrix, 2);
    nLbls = size(lblVec, 2);
    
    areaMatrix = nan(nLbls, nFrames);
    for iFrame = 1:nFrames
       for i=1:nLbls
           lbl = lblVec(i);
           if isempty(cellMatrix{lbl, iFrame})
               curArea = 0;
           else
                curArea = cellMatrix{lbl, iFrame}.area;
           end
           areaMatrix(i, iFrame) = curArea;
       end 
    end
    y = areaMatrix;
end

function y = plotCellPerimeter(cellMatrix, lblVec)
    perimeterFigure = figure;  % this figure is the current active automatically
    perimeterMatrix = getCellPerimeterMatrix(cellMatrix, lblVec);
    nFrames = size(perimeterMatrix, 2);
    plot(1:nFrames, perimeterMatrix);
    title('Perimeter by Frame');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
end

function y = getCellPerimeterMatrix(cellMatrix, lblVec)
    assert(max(lblVec) <= size(cellMatrix, 1));     % make sure all the lbls from the usee are valid
    nFrames = size(cellMatrix, 2);
    nLbls = size(lblVec, 2);
    
    perimeterMatrix = nan(nLbls, nFrames);
    for iFrame = 1:nFrames
       for i=1:nLbls
           lbl = lblVec(i);
           if isempty(cellMatrix{lbl, iFrame})
               curPerimeter = 0;
           else
                curPerimeter = cellMatrix{lbl, iFrame}.perimeter;
           end
           perimeterMatrix(i, iFrame) = curPerimeter;
       end 
    end
    y = perimeterMatrix;
end

% there are 2 types of position- 1 from the original movie and one with
% correction of drift
% todo- i put here 3 types.. need to remove 1
function y = plotCellPosition(cellMatrix, sharedData, lblVec)
    positionFigure = figure;  % this figure is the current active automatically
    %positionMatrix = getCellPositionMatrix(cellMatrix, sharedData, lblVec);
    %nFrames = size(positionMatrix, 2);
    %plot(positionMatrix(:, :, 1).', positionMatrix(:, :, 2).');
    %title('Position by Frame');
    %lblVecLegend = strread(num2str(lblVec),'%s');
    %legend(lblVecLegend);
    %set(gca, 'YDir', 'reverse');
    
    subplot(4,1,1);
    positionMatrix = getCellPositionMatrix(cellMatrix, sharedData, lblVec, 1);
    a = positionMatrix;
    nFrames = size(positionMatrix, 2);
    plot(positionMatrix(:, :, 1).', positionMatrix(:, :, 2).');
    title('Position by Frame- without drift correction');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
    % reverse the aix for the position be align with the movie cordinates
    set(gca, 'YDir', 'reverse');
    
    subplot(4,1,2);
    positionMatrix = getCellPositionMatrix(cellMatrix, sharedData, lblVec, 2);
    b = positionMatrix;
    nFrames = size(positionMatrix, 2);
    plot(positionMatrix(:, :, 1).', positionMatrix(:, :, 2).');
    title('Position by Frame- with drift correction minus');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
    % reverse the aix for the position be align with the movie cordinates
    set(gca, 'YDir', 'reverse');
    
    subplot(4,1,3);
    positionMatrix = getCellPositionMatrix(cellMatrix, sharedData, lblVec, 3);
    b = positionMatrix;
    nFrames = size(positionMatrix, 2);
    plot(positionMatrix(:, :, 1).', positionMatrix(:, :, 2).');
    title('Position by Frame- with drift correction plus');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
    % reverse the aix for the position be align with the movie cordinates
    set(gca, 'YDir', 'reverse');
end

function y = getCellPositionMatrix(cellMatrix, sharedData, lblVec, ver)
    assert(max(lblVec) <= size(cellMatrix, 1));     % make sure all the lbls from the user are valid
    if nargin == 3
        ver = 2;
    end
    nFrames = size(cellMatrix, 2);
    nLbls = size(lblVec, 2);
    
    cordWithDriftCorrectionVer1Minus = nan(nLbls, nFrames, 2);
    cordWithDriftCorrectionVer2Plus = nan(nLbls, nFrames, 2);
    cordNoDriftCorrection = nan(nLbls, nFrames, 2);
    for iFrame = 1:nFrames
       for i=1:nLbls
           lbl = lblVec(i);
           if isempty(cellMatrix{lbl, iFrame})
               currCords = [0 ; 0];
               continue;
           else
               currCords = [cellMatrix{lbl, iFrame}.centroid(1) ; cellMatrix{lbl, iFrame}.centroid(2)];
           end
           a = sum(sharedData.driftPerFrameVer1(1, 1:iFrame));
           b = currCords  - sum(sharedData.driftPerFrameVer1(1, 1:iFrame));
           c = currCords;
           cordWithDriftCorrectionVer1Minus(i, iFrame, :) = currCords  - sum(sharedData.driftPerFrameVer1(1, 1:iFrame));
           cordWithDriftCorrectionVer2Plus(i, iFrame, :) = currCords  + sum(sharedData.driftPerFrameVer1(1, 1:iFrame));
           cordNoDriftCorrection(i, iFrame, :) = currCords;
       end
    end

    if ver==1
        y = cordNoDriftCorrection;
    elseif ver==2
        y = cordWithDriftCorrectionVer1Minus;
    elseif ver==3
        y = cordWithDriftCorrectionVer2Plus;
    elseif ver==4 
    end
end

function y = plotCellNeighbourCounter(cellMatrix, lblVec)
    neighbourCounterFigure = figure;  % this figure is the current active automatically
    neighbourCounterMatrix = getCellNeighbourCounterMatrix(cellMatrix, lblVec);
    nFrames = size(neighbourCounterMatrix, 2);
    plot(1:nFrames, neighbourCounterMatrix);
    ylim([0 10]);
    title('Number of Cell Neighbours by Frame');
    lblVecLegend = strread(num2str(lblVec),'%s');
    legend(lblVecLegend);
end

function y = getCellNeighbourCounterMatrix(cellMatrix, lblVec)
    assert(max(lblVec) <= size(cellMatrix, 1));     % make sure all the lbls from the usee are valid
    nFrames = size(cellMatrix, 2);
    nLbls = size(lblVec, 2);
    
    neighbourCounterMatrix = nan(nLbls, nFrames);
    for iFrame = 1:nFrames
       for i=1:nLbls
           lbl = lblVec(i);
           if isempty(cellMatrix{lbl, iFrame})
               curNeighbour = 0;
           else
                curNeighbour = getNeighborsSum(cellMatrix{lbl, iFrame});
           end
           neighbourCounterMatrix(i, iFrame) = curNeighbour;
       end 
    end
    y = neighbourCounterMatrix;
end

function y = plotCellAreaPerimeterNeighbor(cellMatrix, lblVec)
    areaPerimNeighborFigure = figure;  % this figure is the current active automatically
    lblVecLegend = strread(num2str(lblVec),'%s');
    
    subplot(3,1,1);
    areaMatrix = getCellAreaMatrix(cellMatrix, lblVec);
    nFrames = size(areaMatrix, 2);
    plot(1:nFrames, areaMatrix);
    title('Area by Frame');
    legend(lblVecLegend);
    
    subplot(3,1,2);
    perimeterMatrix = getCellPerimeterMatrix(cellMatrix, lblVec);
    nFrames = size(perimeterMatrix, 2);
    plot(1:nFrames, perimeterMatrix);
    title('Perimeter by Frame');
    legend(lblVecLegend);
    
    subplot(3,1,3);
    neighbourCounterMatrix = getCellNeighbourCounterMatrix(cellMatrix, lblVec);
    nFrames = size(neighbourCounterMatrix, 2);
    plot(1:nFrames, neighbourCounterMatrix);
    title('Number of Cell Neighbours by Frame');
    legend(lblVecLegend);
end

function y = plotSingleCellAreaPerimeterPositionNeighbor(cellMatrix, sharedData, lbl)
    singleCellFigure = figure;  % this figure is the current active automatically

    subplot(2,2,1);
    areaMatrix = getCellAreaMatrix(cellMatrix, lbl);
    nFrames = size(areaMatrix, 2);
    plot(1:nFrames, areaMatrix);
    title('Area by Frame');
    
    subplot(2,2,2);
    perimeterMatrix = getCellPerimeterMatrix(cellMatrix, lbl);
    nFrames = size(perimeterMatrix, 2);
    plot(1:nFrames, perimeterMatrix);
    title('Perimeter by Frame');
    
    subplot(2,2,3);
    neighbourCounterMatrix = getCellNeighbourCounterMatrix(cellMatrix, lbl);
    nFrames = size(neighbourCounterMatrix, 2);
    plot(1:nFrames, neighbourCounterMatrix);
    ylim([0 10]);
    title('Number of Cell Neighbours by Frame');
    
    subplot(2,2,4);
    positionMatrix = getCellPositionMatrix(cellMatrix, sharedData, lbl);
    nFrames = size(positionMatrix, 2);
    plot(positionMatrix(:, :, 1).', positionMatrix(:, :, 2).');
    title('Position by Frame');
    % reverse the aix for the position be align with the movie cordinates
    set(gca, 'YDir', 'reverse');
end

function y = plotAdjacentCellsEdgeSize(cellMatrix, lblVec)
    adjacentCellsEdgeSizeFigure = figure;  % this figure is the current active automatically
    edgeSizeVec = getAdjacentCellsEdge(cellMatrix, lblVec);
    nFrames = size(edgeSizeVec, 2);
    plot(1:nFrames, edgeSizeVec);
    title('Edge Size of 2 Adjacent Cells');
end

function y = getAdjacentCellsEdge(cellMatrix, lblVec)
    assert(max(lblVec) <= size(cellMatrix, 1));     % make sure all the lbls from the usee are valid
    nFrames = size(cellMatrix, 2);
    
    edgeSizeVec = zeros(1, nFrames);
    lbl1 = lblVec(1);
    lbl2 = lblVec(2);
    for iFrame = 1:nFrames
       if isempty(cellMatrix{lbl1, iFrame})
           curCellsEdgeSize = 0;
       else
            curCellsEdgeSize = getNeighborsEdgeSize(cellMatrix{lbl1, iFrame}, lbl2);
       end
       edgeSizeVec(iFrame) = curCellsEdgeSize;
    end
    y = edgeSizeVec;
end

function y = lookForT1Transitions(cellMatrix, handles)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    shortOutFile = strcat('t1-trans-', handles.trackFile);
    outFile = strcat(handles.baseFolder, shortOutFile);
    constrainColorMap = {};
    for iFrame = 2:size(cellMatrix, 2)-4
        for lbl = 1:size(cellMatrix, 1)
            curCell = cellMatrix{lbl, iFrame};
            
            if (isempty(curCell))
                continue;
            end
            if (curCell.isJustCreated == 1)
                continue;
            end
            constrainColorMap = checkCellInT1Transition(cellMatrix, constrainColorMap, curCell);
        end
    end
    colorStackWithConstrains(inFile, outFile, constrainColorMap)
    
end

function y = checkCellInT1Transition(cellMatrix, constrainColorMap, curCell)
    if (curCell.lbl==8 || curCell.lbl==9) && (33 < curCell.frame )
            asd=0;
    end
    
    %disp(strcat('frame: #', num2str(curCell.frame), ', cell #',
    %num2str(curCell.lbl))); todo
    prevCell = cellMatrix{curCell.lbl, curCell.frame-1};
    % we will always have prev because when checked if this cell is not
    % "isJustCreated"
    diffNeighbors1 = setdiff(getNeighborsList(curCell), getNeighborsList(prevCell));
    diffNeighbors1(diffNeighbors1==0) = [];
    if curCell.frame==15 && curCell.lbl==198
        asd=0;
    end
    for diffLbl = diffNeighbors1
        if isT1TrasNeighborAdded(cellMatrix, curCell, diffLbl)
            constrainColorMap = addConstrainColorMap(constrainColorMap, curCell, diffLbl);
        end
    end
    
    diffNeighbors2 = setdiff(getNeighborsList(prevCell), getNeighborsList(curCell));
    diffNeighbors2(diffNeighbors2==0) = [];
    for diffLbl = diffNeighbors2
        if isT1TransNeighborRemoved(cellMatrix, curCell, diffLbl)
            constrainColorMap = addConstrainColorMap(constrainColorMap, curCell, diffLbl);
        end
    end
    y = constrainColorMap;
end

function y = isT1TrasNeighborAdded(cellMatrix, currCell, diffLbl)
    y=1;
    if currCell.frame==8 && currCell.lbl==211
        asd=0;
    end
    if cellMatrix{diffLbl, currCell.frame}.isJustCreated == 1
        y = 0;
    else
        for iFrame = currCell.frame+1:currCell.frame+3
            c = cellMatrix{currCell.lbl, iFrame};
            if isempty(c)
                y = 0;
                break;
            end
            if not(any(getNeighborsList(c)==diffLbl))
                %  if this Neighbor disppeare from list (after it was
                %  added)- this is not t1
                y = 0;
                break;
            end
        end
    end
end

function y = isT1TransNeighborRemoved(cellMatrix, currCell, diffLbl)
    y=1;
    if isempty(cellMatrix{diffLbl, currCell.frame})
        % if the diff is because a cell that  just died
        y = 0;
    else
        for iFrame = currCell.frame+1:currCell.frame+3
            c = cellMatrix{currCell.lbl, iFrame};
            if isempty(c)
                y = 0;
                break;
            end
            if isempty(cellMatrix{diffLbl, iFrame})
                y=0;
                break;
            end
            if any(getNeighborsList(c)==diffLbl)
                y = 0;
                break;
            end
        end
    end
end

% this function create a new movie- with constains- each constrain is a
% 3-tupele contain- frame no. and 2 lable, and it should print these lbl in
% the frame in red
function y = addConstrainColorMap(constrainColorMap, curCell, diffLbl)
    iConstrain = size(constrainColorMap, 1) + 1;
    constrainColorMap{iConstrain, 1} = curCell.frame;
    constrainColorMap{iConstrain, 2} = [curCell.lbl diffLbl];
    y = constrainColorMap;
end

% create new movie from constrain map- please see addConstrainColorMap()
function y = colorStackWithConstrains(inFile, outFile, constrainColorMap)
    if exist(outFile, 'file')
        delete(outFile);
    end
    
    n = 500;
    randVec = rand(n, 1);
    baseColorMap = [randVec randVec randVec];
    
    stackInfo = imfinfo(inFile);
    nFrames = numel(stackInfo);

    for iFrame = 1:nFrames
      curFrame = imread(inFile, iFrame, 'Info', stackInfo);
      colorMap = buildColorMap(baseColorMap, constrainColorMap, iFrame);
      ColoredFrame = label2rgb(curFrame, colorMap , 'k');
      imwrite(uint16(ColoredFrame), outFile, 'WriteMode', 'append',  'Compression','none');
    end
end

function y = buildColorMap(baseColorMap, constrainColorMap, iFrame)
    for i = 1:size(constrainColorMap, 1)
        if iFrame-1 <= constrainColorMap{i,1} && constrainColorMap{i,1} <= iFrame+1
            curColor = i * 0.08;
            for lbl = constrainColorMap{i,2}
                %baseColorMap(lbl, :) = [curColor 0 0];
                baseColorMap(lbl, :) = [1 0 0];
            end
        end
    end
    y = baseColorMap;
end

function y = plotCellOnStack(lblVec, handles)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    shortOutFile = strcat('single-color-', handles.trackFile);
    outFile = strcat(handles.baseFolder, shortOutFile);
    n = 500;
    randVec = rand(n, 1);
    myColorMap = [randVec randVec randVec]; % creare colorMap of grayscale
    for lbl=lblVec
        myColorMap(lbl, :) = [1 0 0];
    end
    
    colorStack(inFile, outFile, myColorMap)
end

function y = colorStackRandom(handles)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    shortOutFile = strcat('color-', handles.trackFile);
    outFile = strcat(handles.baseFolder, shortOutFile);
    n = 500;
    myColorMap = rand(n, 3);
    colorStack(inFile, outFile, myColorMap)
end

% create new movie with give colormap
function y = colorStack(inFile, outFile, colorMap)
    initFile(outFile);
    stackInfo = imfinfo(inFile);
    nFrames = numel(stackInfo);

    for iFrame = 1:nFrames
      curFrame = imread(inFile, iFrame, 'Info', stackInfo);
      ColoredFrame = label2rgb(curFrame, colorMap , 'k');
      imwrite(uint16(ColoredFrame), outFile, 'WriteMode', 'append',  'Compression','none');
    end
end

function y = menualCorrection(handles, lblVec)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    shortOutFile = strcat('menual-correction-', handles.trackFile);
    outFile = strcat(handles.baseFolder, shortOutFile);
    
    initFile(outFile);
    stackInfo = imfinfo(inFile);
    nFrames = numel(stackInfo);
    newLbl = lblVec(1);
    for iFrame = 1:nFrames
       oldFrame = imread(inFile, iFrame, 'Info', stackInfo);
       fixedFrame = oldFrame;
       for lbl = lblVec(2:end)
           pixelOfLbl = find(oldFrame==lbl);
           fixedFrame(pixelOfLbl) = newLbl;
       end
       imwrite(uint16(fixedFrame), outFile, 'WriteMode', 'append',  'Compression','none');
    end
end

% interactive choose lbl for other cmd
% this will listen for the user and generate lbl vector for the other cmds
function y = chooseLbl(hObject, handles, iFrame)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    stackInfo = imfinfo(inFile);
    curFrame = imread(inFile, iFrame, 'Info', stackInfo);
    chooseLblFigure = figure;  % this figure is the current active automatically
    imshow(curFrame, []);
    [yList, xList] = getpts;    % listen to the user mouse
    close(chooseLblFigure);
    lblVec = [];
    for i = 1:size(xList, 1)
       x = uint32(xList(i));
       y = uint32(yList(i));
       newLbl = curFrame(x,y);
       lblVec = [lblVec curFrame(x,y)];
    end
    UniqueIndx = unique(lblVec);
    UniqueIndx(UniqueIndx==0) = [];
    lblVec = strjoin(string(UniqueIndx));
    set(handles.text_input, 'String', lblVec);  % insert the clicked cells into the GUI textbox
end

% generate a new movie that fix the global drift of the cells
function fixDrift(handles)
    inFile = strcat(handles.baseFolder, handles.trackFile);
    shortOutFile = strcat('fixe-drift-', handles.trackFile);
    outFile = strcat(handles.baseFolder, shortOutFile);
    
    initFile(outFile);
    stackInfo = imfinfo(inFile);
    nFrames = numel(stackInfo);
    for iFrame = 1:nFrames
       curFrame = imread(inFile, iFrame, 'Info', stackInfo);
       rowOffset = uint32(sum(handles.sharedData.driftPerFrameVer1(1, 1:iFrame)))*1;
       colOffset = uint32(sum(handles.sharedData.driftPerFrameVer1(2, 1:iFrame)))*1;
       
       fixedFrame = pasteOver(curFrame, double(colOffset), double(rowOffset));
       imwrite(uint16(fixedFrame), outFile, 'WriteMode', 'append',  'Compression','none');
    end
    disp('Done');
end

% colOffset, rowOffset are the direction the movie went- we will offset the
% new movie to the other direction
function y = pasteOver(small, colOffset, rowOffset)
    [rowSize, colSize] = size(small);
    constOffset = 150;
    padBefore.row = constOffset+rowOffset;
    padBefore.col = constOffset+colOffset;
    
    padAfter.row = constOffset-rowOffset;
    padAfter.col = constOffset-colOffset;
    
    newImage = small;
    newImage = padarray(newImage, [padBefore.row padBefore.col], 'pre');
    newImage = padarray(newImage, [padAfter.row padAfter.col], 'post');
    y = newImage;
end

function y = initFile(fileName)
    if exist(fileName, 'file')
        delete(fileName);
    end
end