%clear all;
%close all;

main();

function y = main()
    disp('Step #1 Started');
    WinAreaFrc = 0.7;       % Predominant area fraction
    ShareAreaFrc = 0.4;     % Significant area fraction
    eyalShareAreaFrc = 0.55;
    minCellSize = 52;
    backTraceFactor = 0.2;
    backTraceMinFactor = 0.5;
    posibleSplitRatio = 0.8;

    [movieFile,baseFolder] = uigetfile('*.tif','select the segmented movie to analyse');
    inMovieFile = strcat(baseFolder, movieFile);
    if not(exist(inMovieFile, 'file'))
        uiwait(msgbox('input movie was not selected- cmd aborted', 'Error','error'));
        return;
    end

    testFile = strcat(baseFolder, 'test.tif');
    if exist(testFile, 'file')
        uiwait(msgbox('testMovie was not deleted', 'Error','error'));
        return;
    end
    
    outFileName = strcat(baseFolder,strcat('out-', movieFile));
    initFile(outFileName);

    stackInfo = imfinfo(inMovieFile);
    nFrames = numel(stackInfo);

    for iFrame = 1:nFrames
        curFrame = imread(inMovieFile, iFrame, 'Info', stackInfo);
        [imageHeight, imageWidth] = size(curFrame);
        curFrame = imclearborder(255-curFrame);

        CC = bwconncomp(curFrame==255, 8);
        Particles = CC.PixelIdxList;    
        nPart = length(Particles);

        curTmpLblMap = zeros(size(curFrame));
        AreaPart = zeros(1,nPart);
        for i = 1:nPart
            curTmpLblMap(Particles{i}) = i; 
            AreaPart(i) = length(Particles{i});
        end

        if iFrame==1
            lblMap = curTmpLblMap;
            lblCounter = nPart;

            isWeMergeCell = zeros(1, nPart*4);
        else
            imwrite(uint16(curFrame), testFile, 'WriteMode', 'append',  'Compression','none');
            % prevLblMap will help us figure for the current cell which lbl
            % we should give in- according to the prev image
            prevLblMap = lblMap.*(curTmpLblMap>0);

            % for each current cell- this cell vector store all the lbls
            % fro mprev frame that are on the same area as this cell
            cellLblAndArea = cell(1,nPart);
            for i = 1:nPart
                Labels = prevLblMap(Particles{i}).';
                UniqueIndx = unique(Labels);
                UniqueIndx(UniqueIndx==0) = [];
                AreaVector = zeros(1,length(UniqueIndx));
                for k =1:length(UniqueIndx)
                    msk = (Labels==UniqueIndx(k));
                    AreaVector(k) = sum(msk);
                end
                cellLblAndArea{i} = sortrows([UniqueIndx ; AreaVector]',2, 'descend')';
            end

            lblMap = prevLblMap;
            curLbl2MainLbl = zeros(nPart, 3);

            %% assign just lbl we are sure of
            for i = 1:nPart
                if(not(isempty(cellLblAndArea{i})))
                    [mxIndx, mxIndxArea] = getMaxIndexAndArea(cellLblAndArea{i});

                    % is the cell has been split to several cells-
                    % e.g. border was added- dont handle it now
                    taegetLblPrevArea = nnz(prevLblMap==mxIndx);
                    cond1 = mxIndxArea/taegetLblPrevArea > posibleSplitRatio;

                    % if most of this cell area contains the mxIndx lbl-
                    % this lbl is probebly a good match
                    cond2 = mxIndxArea/AreaPart(i)>WinAreaFrc;
                    cond3 = mxIndxArea>ShareAreaFrc*AreaPart(i) && length(cellLblAndArea{i}(1,:))==1;

                    if (cond1 && (cond2 || cond3))
                        lbl = mxIndx;
                        assignPixel = Particles{i};
                        % assign lbl
                        lblMap(assignPixel) = lbl;

                        tmp = find(curLbl2MainLbl(i,:)==0, 1);
                        curLbl2MainLbl(i, tmp) = lbl;
                        isWeMergeCell(lbl) = isWeMergeCell(lbl)*backTraceMinFactor;
                        % end assign lbl
                    end
                end   
            end

            %% border disapeared
            % check if need to split the current cell into 2- because in the
            % last frame this cell was 1 cells
            for i = 1:nPart
                if (not(curLbl2MainLbl(i, 1) == 0))
                    continue;
                end
                if (getNumOfLblOnCell(cellLblAndArea{i}) < 2)
                    % if there are less that 2 lbl in cell- not split can
                    % occure...
                    continue;
                end
                if (sum(AreaPart(i)) < minCellSize)
                    % if cell is too small- 
                    % probably a border was added and not disapeared..
                    continue;
                end
                if (sum(cellLblAndArea{i}(2,:))*1.2 < AreaPart(i))
                    % if most of this cell area is new- we probebly need to
                    % make this a brand new cell, and now split it
                    continue;
                end
    
                % check if we can reasonably split this cell into 2/3
                if isCellSplitIntoX(cellLblAndArea{i}, 3)
                    lblSplitted = cellLblAndArea{i}(1, 1:3);
                elseif isCellSplitIntoX(cellLblAndArea{i}, 2)
                    lblSplitted = cellLblAndArea{i}(1, 1:2);
                else
                    % this cell cannot be splitted
                    continue;
                end

                allHandlesSplittedCells = [];   % this is a list of all this cell pixel that we assign a lbl
                centroVec = [];                 % this is a vector of all the centroid of the lbl this cell will get- we use this to calc the lbl for pixel that are not on any of the prev lbls
                for j = 1:size(lblSplitted, 2)
                    lbl = lblSplitted(j);
                    originalCellPixels = Particles{i};
                    prevLblPixels = find(prevLblMap==lbl);
                    assignPixel = intersect(originalCellPixels, prevLblPixels);
                    % assign lbl
                    lblMap(assignPixel) = lbl;

                    tmp = find(curLbl2MainLbl(i,:)==0, 1);
                    curLbl2MainLbl(i, tmp) = lbl;
                    isWeMergeCell(lbl) = isWeMergeCell(lbl)*backTraceMinFactor;
                    % end assign lbl

                    allHandlesSplittedCells = vertcat(allHandlesSplittedCells, assignPixel);

                    % next we remove small borders of this cell, so we will
                    % have just 1 centroid for each cell
                    seSplit = strel('disk',10);
                    splitCentroidImg = imclose(prevLblMap==lbl, seSplit);
                    stats1 = regionprops(splitCentroidImg,'Centroid');
                    
                    curCentroid = stats1.Centroid;
                    centroVec = vertcat(centroVec, curCentroid);
                end
 
                % we will have lbls that are not in any of the prev lbls-
                % so we will calc their closest cell, and give then this
                % lbl
                unhandledIndex = setdiff(Particles{i}, allHandlesSplittedCells);
                for j = 1:length(unhandledIndex)
                    curPixel = unhandledIndex(j);
                    curCord = [curPixel/imageHeight, mod(curPixel,imageHeight)];
                    closestCellIndex = calcClosestCell(curCord, centroVec);

                    lbl = lblSplitted(closestCellIndex);
                    assignPixel = curPixel;
                    % assign lbl
                    lblMap(assignPixel) = lbl;
                    % end assign lbl
                end
            end


            %% make or abort split (border was added)
            % generally we abort split- but uf we see this border is showing
            % again- we will make the split

            % in the prev frame this cell was part of a larger cell-
            % here we check if this large cell is segmentation error
            % or is it a permanent state- so we should create a seperate cell
            % for each new lbl
            for i = 1:nPart
                if (not(curLbl2MainLbl(i, 1) == 0))
                    continue;
                end
                if (isempty(cellLblAndArea{i}))
                    continue;
                end

                [mxIndx, mxIndxArea] = getMaxIndexAndArea(cellLblAndArea{i});
                if (1 < isWeMergeCell(mxIndx))
                    % apperently this cell has turned into 1- segmentation error
                    % we need to abort the merge and make this a seperate single cell
                    % for each lbl that compound the prev large cell we create
                    % new seperate cell
                    lblCounter = lblCounter+1;
                    lbl = lblCounter;
                    assignPixel = Particles{i};

                    % assign lbl
                    lblMap(assignPixel) = lbl;

                    tmp = find(curLbl2MainLbl(i,:)==0, 1);
                    curLbl2MainLbl(i, tmp) = lbl;
                    isWeMergeCell(lbl) = 0;
                    % end assign lbl
                end
            end

            for i = 1:nPart
                if (not(curLbl2MainLbl(i, 1) == 0))
                    continue;
                end
                if (isempty(cellLblAndArea{i}))
                    continue;
                end

                [mxIndx, mxIndxArea] = getMaxIndexAndArea(cellLblAndArea{i});
                if (eyalShareAreaFrc < mxIndxArea/AreaPart(i))
                    % we assume the new border was added by accedent so we
                    % eliminate it and merge the cell into a bigger one..
                    % (because in the previos frames this was the case)
                    % in addition, in case we are wrong we will increase
                    % the isWeMergeCell for this main lbl
                    lblMap(Particles{i}) = mxIndx;

                    tmp = find(curLbl2MainLbl(i,:)==0, 1);
                    curLbl2MainLbl(i, tmp) = mxIndx;
                    isWeMergeCell(mxIndx) = isWeMergeCell(mxIndx) + backTraceFactor;
                end
            end


            % Find appearing object (particles not assigned any object label)
            % Avoid artefacts during division
            se2=strel('disk',2);
            ObjLblDilate = imdilate(lblMap,se2);
            DiffMask = ((curTmpLblMap>0)-(ObjLblDilate>0))>0;
            if any(any(DiffMask,1))
                CC2 = bwconncomp(DiffMask, 8);
                for i=1:CC2.NumObjects
                    cellSize = size(CC2.PixelIdxList{i}, 1);
                    if cellSize < minCellSize
                        %continue;
                    end
                    % create new lbl and assign this lbl to this new cell
                    lblCounter = lblCounter+1;
                    lbl = lblCounter;
                    curLbl = unique(curTmpLblMap(CC2.PixelIdxList{i}));
                    curLbl(curLbl==0) = [];
                    assignPixel = Particles{curLbl};

                    % assign lbl
                    lblMap(assignPixel) = lbl;

                    tmp = find(curLbl2MainLbl(curLbl,:)==0, 1);
                    curLbl2MainLbl(curLbl, tmp) = lbl;
                    isWeMergeCell(lbl) = 0;
                    % end assign lbl
                end
            end

            %% Merge touching subparticles assigned to the same object label
            se=strel('disk',1);
            ObjLblClosed = imclose(lblMap,se);
            for i=1:lblCounter
                CC = bwconncomp(ObjLblClosed==i, 8);
                for j=1:CC.NumObjects
                    lblMap(CC.PixelIdxList{j}) = i;
                end
            end

            %% Redraw object boundaries
            se=strel('disk',1);
            ObjDilate = imdilate(lblMap,se);
            ObjLblBnd = lblMap;
            ObjLblBnd(ObjDilate~=lblMap) = 0;
            lblMap = ObjLblBnd;

            imwrite(uint16(lblMap), outFileName, 'WriteMode', 'append',  'Compression','none');
        end  
    end
    disp('Step #1 Ended');
    uiwait(msgbox(strcat('Step #1 Done- movie was analysed and cell were tracked, output file is ', outFileName)));
end

function y = initFile(fileName)
    if exist(fileName, 'file')
        delete(fileName);
    end
end

function [indx area] = getMaxIndexAndArea(cellLblAndArea)
    indx = cellLblAndArea(1,1);
    area= cellLblAndArea(2,1);
end

function y = getNumOfLblOnCell(cellLblAndArea)
    y = size(cellLblAndArea, 2);    % return 2nd dimention size
end


function y = isCellSplitIntoX(cellLblAndArea, x)
    if (size(cellLblAndArea) < x)
        y = false;
    else
        MIN_CELL_SIZE = 100;
        eyalFrc = 0.82;
        % count how many area are more that 50 pixels
        areaVec = cellLblAndArea(2,:);
        y = sum(areaVec>MIN_CELL_SIZE);

        result1 = (y == x);
        areaRatio = sum(areaVec(1:x)) / sum(areaVec);
        result2 = areaRatio > eyalFrc;
        y = result1 && result2;
    end
end


function y = calcClosestCell(curCord, centroVec);
    closestIndex = 1;
    closestVal = calcDistance(curCord, centroVec(1,:));
    for i = 2:size(centroVec, 1)
        curDistance = calcDistance(curCord, centroVec(i,:));
        if (curDistance < closestVal)
            closestIndex = i;
            closestVal = curDistance;
        end
    end
    y = closestIndex;
end

function y = calcDistance(cordA, cordB)
    x1 = cordA(1);
    y1 = cordA(2);
    x2 = cordB(1);
    y2 = cordB(2);
    
    % no need to make the sqrt since this is symetric to akk distance- it
    % will save us time
    y = (x2-x1).^2 + (y2-y1).^2;
end

function y = isNeedToAbortMerge(cellLblAndArea, isWeMergeCell)
    
    if (size(cellLblAndArea, 2) == 0)
        eyalAssert();
    elseif (size(cellLblAndArea, 2) == 1)
        y = isWeMergeCell(cellLblAndArea(1));
    else
        y = isWeMergeCell(cellLblAndArea(1));
        for i = 2:size(cellLblAndArea, 2)
            x = isWeMergeCell(cellLblAndArea(1));
            if (y == true && not(x==y))
                eyalassert();
            end
        end
    end
end

function y = eyalAssert()
    eyalAssertTmp = 4;
end