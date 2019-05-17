%clear all;
%close all;

% this is the main function for step #2
function y = main()
    disp('Step #2 Started');
    [movieFile,baseFolder] = uigetfile('*.tif','select the segmented movie to parse for DB');
    stackFile = strcat(baseFolder, movieFile);
    if not(exist(stackFile, 'file'))
        uiwait(msgbox('input movie was not selected- abort cmd', 'Error','error'));
        return;
    end
    
    csvFile = strcat(baseFolder, 'csvCell.csv');

    cellMatrix = fillCellMatrix(stackFile); % from the segmented movie create the DB
    writeCellMatrixToCsv(cellMatrix, csvFile); % store the DB to a file
    disp('Step #2 Finished');
    uiwait(msgbox('Step #2 Done- Movie DB was created. the DB file is stored in ./csvCell.csv'));
end

% for each cell in each frame create a new MyCell instance and push it to
% the cell matrix
function y = fillCellMatrix(stackFile)
    stackInfo = imfinfo(stackFile);
    nFrames = numel(stackInfo);
    
    lastFrame = imread(stackFile, nFrames, 'Info', stackInfo);
    lablesInLastFrame = unique(lastFrame);
    maxLbl2 = max(lablesInLastFrame);
    cellMatrix = cell(maxLbl2, nFrames);
    for iFrame = 1:nFrames
       curFrame = imread(stackFile, iFrame, 'Info', stackInfo);
       lablesInFrame = unique(curFrame);
       lablesInFrame(lablesInFrame==0) = [];
       
       for i = 1:size(lablesInFrame,1)
           lbl = lablesInFrame(i);
           c = MyCell;
           c.frame = iFrame;
           c.lbl = lbl;           
           Obj = (curFrame == lbl);
           
           c.area = sum(sum(Obj));
           statsStruct = regionprops(Obj,{'Centroid', 'perimeter'});
           c.perimeter = statsStruct.Perimeter;
           
           c.centroid = statsStruct.Centroid;
           
           c.isJustCreated = 1;
           c.prevCentroid = c.centroid;
           if ((1 < c.frame) && (c.lbl <= size(cellMatrix, 1)))
               if not(isempty(cellMatrix{c.lbl, c.frame-1}))
                    c.prevCentroid = cellMatrix{c.lbl, c.frame-1}.centroid;
                    c.isJustCreated = 0;
                    % if we get here- then this cell has predecessor and
                    % this is not the first frame..
                    % so this cell was NOT just created
               end
           end
           c.driftDistance = -1;   % this will be calc later on 
           c.driftAngle = -1;      % this will be calc later on
           c.neighborsList = getNeighbours(curFrame, c);
           
           cellMatrix{c.lbl, c.frame} = c;
       end
       fprintf('frame #%d finished', iFrame);
    end
    assert(size(cellMatrix, 2) == nFrames);
    
    y = cellMatrix;
end

% create .csv (Comma-separated values) file, each row represent a cell in
% some frame
function y = writeCellMatrixToCsv(cellMatrix, csvFile)
    initFile(csvFile);
    for iFrame = 1:size(cellMatrix, 2)
        for lbl = 1:size(cellMatrix, 1)
            curCell = cellMatrix{lbl, iFrame};
            if (isempty(curCell))
                continue;
            end
            dlmwrite(csvFile, toMatrix(curCell), '-append');
        end
        disp(strcat("frame #", num2str(iFrame), " finished"));
    end
end

% use morphological operations to generate neighbours list for each cell
function y = getNeighbours(curFrame, curCell)
    object = (curFrame == curCell.lbl);
    se = ones(8);
    neighbours = imdilate(object, se) & ~object;
    neighbourLabels = unique(curFrame(neighbours));
    neighbourLabels(neighbourLabels==0) = [];
    neighbourLabels = neighbourLabels.';
    
    neighboursAndBorders = [];
    if size(neighbourLabels, 1) == 0
        neighboursAndBorders=[0;0];
    else
        for lbl = neighbourLabels
            currNeighbourAndBorder = [lbl ; getBorderSize(curFrame, curCell.lbl, lbl, se)];
            neighboursAndBorders = [neighboursAndBorders currNeighbourAndBorder];
        end
    end
    
    % the structure of the .csv file is defined well so we need to make
    % sure the neighbours list is exactly in size 10
    if size(neighboursAndBorders, 2) < 10
        elementsToPad = 10-size(neighboursAndBorders, 2);
        neighboursAndBorders = padarray(neighboursAndBorders, [0 elementsToPad], 'post');
    else
        neighboursAndBorders = neighboursAndBorders(:, 1:10)
    end
    
    y = neighboursAndBorders;
end

function y = getBorderSize(curFrame, lbl1, lbl2, se)
    object1 = (curFrame == lbl1);
    lblBorder1 = imdilate(object1, se) & ~object1;

    object2 = (curFrame == lbl2);
    lblBorder2 = imdilate(object2, se) & ~object2;

    finalBorder = lblBorder1 & lblBorder2;
    y = nnz(finalBorder);
end

function y = initFile(fileName)
    if exist(fileName, 'file')
        delete(fileName);
    end
end