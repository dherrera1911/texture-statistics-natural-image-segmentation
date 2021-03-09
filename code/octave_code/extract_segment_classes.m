%%%%%%% Generate statistics data from texture database %%%%%%%%%
pkg load image
pkg load io
addpath(genpath("./functions"));

imageDir = "../../images/BSR/BSDS500/data/images/all/";
segmentsDir = "../../images/BSR/BSDS500/data/groundTruth/all/"
classificationDir = "../../data/BSD_results/6_classification_agreement.csv"
savingDir = "../../data/sorted_images/";

# read contents of a dir that are not hidden files
function [files] = read_unhidden(inputDir)
  temp = dir(inputDir);
  temp = {temp.name};
  hiddenFiles = regexp(temp, "^\\.");
  unhiddenFiles = cellfun("isempty", hiddenFiles);
  files = temp(unhiddenFiles);
end

% return whether the strings in a cell are an image file type
function [isImage] = is_image_file(strCell)
  imageTypes = {".tif", ".jpg", ".png", ".tiff", ".pnm"};
  isImage = zeros(length(strCell), 1);
  for it = 1:length(imageTypes)
    matchInd = strfind(lala, imageTypes{it});
    matchInd = !(cellfun("isempty", matchInd));
    isImage = isImage + matchInd;
  end
end

%erenurn fsrst value of vector different from reference value
function outVal = first_different(inputVec, referenceVal)
  ind = find(inputVec != referenceVal);
  if isempty(ind)
    outVal = NA;
  else
    outVal = double(inputVec(ind(1)));
  end
end

%{
Get a segmentation map as input, and initial coords for a given pixel.
Return a vector indicating the neighboring segments up, down, left and
right, returning NA if no neighboring segment in that direction.
initialCoords are in [row col]
%}
function neighborSegments = find_neighbor(segmentMap, initialCoords)
  xi = initialCoords(2);
  yi = initialCoords(1);
  height = size(segmentMap)(1);
  width = size(segmentMap)(2);
  mainSeg = segmentMap(yi, xi);
  mainSegMask = segmentMap == mainSeg;
  % find neighbouring segments by applying offsets
  % offset left
  segmentsL = unique(segmentMap(logical([mainSegMask(:,2:end), zeros(height,1)])));
  segmentsR = unique(segmentMap(logical([zeros(height,1), mainSegMask(:,1:(end-1))])));
  segmentsD = unique(segmentMap(logical([zeros(1,width); mainSegMask((1:(end-1)),:)])));
  segmentsU = unique(segmentMap(logical([mainSegMask((2:end),:); zeros(1,width)])));
  neighborSegments = unique([segmentsL', segmentsR', segmentsU', segmentsD']);
  neighborSegments = neighborSegments(neighborSegments != mainSeg);
end

%{
Filter segments to have a minimum size
%}
function filteredSegments = filter_segments(neighborSegments, segmentMap, minSize = 64*64)
  filteredSegments = [];
  for seg = neighborSegments
    nPix = sum(segmentMap(:) == seg);
    if (nPix >= minSize)
      filteredSegments = [filteredSegments, seg];
    end
  end
end 

%{
Find the horizontal centroid of a mask of 0's and 1's
%}
function centroidCoord = find_centroid(inputMask)
  horizontalSum = sum(inputMask);
  hCumulative = cumsum(horizontalSum);
  hCumulative = hCumulative / hCumulative(end) - 0.5;
  centroidCoord = find(hCumulative >= 0)(1);
end

%{
Cut mask in half by horizontal centroic
%}
function cutMask = split_mask(inputMask)
  xHalf = find_centroid(inputMask);
  cutMask{1} = inputMask;
  cutMask{1}(:, xHalf:end) = 0;
  cutMask{2} = inputMask;
  cutMask{2}(:, 1:(xHalf-1)) = 0;
end

%%%%%%%%%%%% Load the classifications for the pairs of segments %%%%%%%%%%%%
classCell = csv2cell(classificationDir);
colNames = classMat(1,:);
designInds = find(ismember(colNames, {"same", "ImageName", "Type", "Segment", ...
"betterHOS", "betterFASHOS"}));
designCell = classCell(:, designInds);
classStruct = cell2struct(designCell(2:end,:).', designCell(1,:));
for i = 1:length(classStruct)
  classStruct(i).ImageName = num2str(classStruct(i).ImageName);
end

%%%%%%%%%%%%
imageFiles = read_unhidden(imageDir);
segmentFiles = read_unhidden(segmentsDir);

minSize = 64*64;

ti = 0;
hi = 0;
ni = 0;
nhi = 0;

for im = 1:length(imageFiles)
  imName = strsplit(imageFiles{im}, "."){1};
  % load image and segmentation map
  imageBSD = imread([imageDir, imageFiles{im}]);
  imageBSD = rgb2gray(imageBSD);
  segment = load([segmentsDir, segmentFiles{im}]);
  segmentMap = segment.groundTruth{2}.Segmentation;
  % remove 1 column and row to make conformant with steerable pyramid
  imageBSD = imageBSD(1:(end-1), 1:(end-1));
  segmentMap = segmentMap(1:(end-1), 1:(end-1));
  %get sizes and center coords
  imageSize = size(imageBSD);
  xCenter = round(imageSize(2)/2);
  yCenter = round(imageSize(1)/2);
  % get center segment and neighbor segments
  centerSeg = segmentMap(yCenter, xCenter);
  neighborSegments = unique(find_neighbor(segmentMap, [yCenter, xCenter]));
  neighborSegments = neighborSegments(!isna(neighborSegments)); 
  % filter by minimum size
  neighborSegments = filter_segments(neighborSegments, segmentMap, minSize);
  centerMask = double(segmentMap == centerSeg);
  sizeCenter = sum(centerMask(:));
  # extract data for this image from class struct
  imageStruct = classStruct(find(ismember({classStruct.ImageName}, {imName})));

  if (!isempty(neighborSegments) & (sizeCenter >= minSize*2))
    ti = ti + 1;
    % center stats
    centerImage{ti} = imName;
    % split center stats
    #splitMask = split_mask(centerMask);
    ##########################
    # Center segment
    ##########################
    # extract center class
    centerStruct = imageStruct(find(ismember({imageStruct.Type}, {"half_center"})));
    # Save in directory corresponding to center class
    if (centerStruct.betterHOS==1)
      centerHOSFile = [savingDir, "betterHOS/", imName, "_center_same.png"];
    else
      centerHOSFile = [savingDir, "nonBetterHOS/", imName, "_center_same.png"];
    end
    if (centerStruct.betterFASHOS==1)
      centerFASHOSFile = [savingDir, "betterFASHOS/", imName, "_center_same.png"];
    else
      centerFASHOSFile = [savingDir, "nonBetterFASHOS/", imName, "_center_same.png"];
    end
    imwrite(imageBSD, centerHOSFile, 'Alpha', 255*centerMask);
    imwrite(imageBSD, centerFASHOSFile, 'Alpha', 255*centerMask);
    ##########################
    # Neighbor segments
    ##########################
    % neighbor stats
    neighborStruct = imageStruct(find(ismember({imageStruct.Type}, {"neighbor"})));
    halfNeighborStruct = imageStruct(find(ismember({imageStruct.Type}, {"half_neighbor"})));
    for n = 1:length(neighborSegments)
      neighborSegment = n;
      segmentStruct = neighborStruct(n);
      neighborMask = double(segmentMap == neighborSegments(n)); 
      pairMask = neighborMask + centerMask;
      % Save in directory corresponding to pair class
      if (segmentStruct.betterHOS==1)
        neighborHOSFile = [savingDir, "betterHOS/", imName, "_n", ...
          num2str(n), "_different.png"];
      else
        neighborHOSFile = [savingDir, "nonBetterHOS/", imName, "_n", ...
          num2str(n), "_different.png"];
      end
      if (segmentStruct.betterFASHOS==1)
        neighborFASHOSFile = [savingDir, "betterFASHOS/", imName, "_n", ...
          num2str(n), "_different.png"];
      else
        neighborFASHOSFile = [savingDir, "nonBetterFASHOS/", imName, "_n", ...
          num2str(n), "_different.png"];
      end
      imwrite(imageBSD, neighborHOSFile, 'Alpha', 255*pairMask);
      imwrite(imageBSD, neighborFASHOSFile, 'Alpha', 255*pairMask);
      % if large enough, compute statistics of half-neighbor
      if (sum(neighborMask(:)) >= minSize*2)
        segmentStruct = halfNeighborStruct(find([halfNeighborStruct.Segment]==n));
        %splitMask = split_mask(neighborMask);
        if (segmentStruct.betterHOS==1)
          halfNeighborHOSFile = [savingDir, "betterHOS/", imName, "_n", ...
            num2str(n), "_neighborSame.png"];
        else
          halfNeighborHOSFile = [savingDir, "nonBetterHOS/", imName, "_n", ...
            num2str(n), "_neighborSame.png"];
        end
        if (segmentStruct.betterFASHOS==1)
          halfNeighborFASHOSFile = [savingDir, "betterFASHOS/", imName, "_n", ...
            num2str(n), "_neighborSame.png"];
        else
          halfNeighborFASHOSFile = [savingDir, "nonBetterFASHOS/", imName, "_n", ...
            num2str(n), "_neighborSame.png"];
        end
        imwrite(imageBSD, halfNeighborHOSFile, 'Alpha', 255*neighborMask);
        imwrite(imageBSD, halfNeighborFASHOSFile, 'Alpha', 255*neighborMask);
      end
    end
  end 
end

