%%%%%%% Generate statistics data from texture database %%%%%%%%%
pkg load image
addpath(genpath("./functions"));

imageDir = "../../images/BSR/BSDS500/data/images/all/";
segmentsDir = "../../images/BSR/BSDS500/data/groundTruth/all/"
savingDir = "../../data/BSD_stats/";

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

%%%%%%%%%%%% Compute the stats %%%%%%%%%%%%%%%%
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

  if (!isempty(neighborSegments) & (sizeCenter >= minSize*2))
    ti = ti + 1;
    % center stats
    centerImage{ti} = imName;
    opts = metamerOpts(zeros(imageSize), 'windowType=flexible', 'verbose=0', 'Na=7');
    opts.flexibleMasks(1,:,:) = centerMask;
    masks = mkFlexibleMasks(opts);
    paramsTemp = metamerAnalysisCorr(double(imageBSD), masks, opts);
    paramsTemp.oim = NA;
    paramsCenter(ti) = paramsTemp; 
    centerImage{ti} = imName;
    % split center stats
    splitMask = split_mask(centerMask);
    for h = 1:2
      hi = hi + 1;
      halfCenterImage{hi} = imName;
      halfCenterNumber{hi} = h;
      % stats of the center segment separated in half
      opts = metamerOpts(zeros(imageSize), 'windowType=flexible', 'verbose=0', 'Na=7');
      opts.flexibleMasks(1,:,:) = splitMask{h};
      masks = mkFlexibleMasks(opts);
      paramsTemp = metamerAnalysisCorr(double(imageBSD), masks, opts);
      paramsTemp.oim = NA;
      paramsCenterHalf(hi) = paramsTemp; 
    end
    % neighbor stats
    for n = 1:length(neighborSegments)
      ni = ni + 1;
      neighborImage{ni} = imName;
      neighborSegment{ni} = n;
      neighborMask = double(segmentMap == neighborSegments(n)); 
      opts = metamerOpts(zeros(imageSize), 'windowType=flexible', 'verbose=0', 'Na=7');
      opts.flexibleMasks(1,:,:) = neighborMask;
      masks = mkFlexibleMasks(opts);
      paramsTemp = metamerAnalysisCorr(double(imageBSD), masks, opts);
      paramsTemp.oim = NA;
      paramsNeighbor(ni) = paramsTemp; 
      % if large enough, compute statistics of half-neighbor
      if (sum(neighborMask(:)) >= minSize*2)
        splitMask = split_mask(neighborMask);
        for h = 1:2
          nhi = nhi + 1;
          halfNeighborImage{nhi} = imName;
          halfNeighborSegment{nhi} = n;
          halfNeighborNumber{nhi} = h;
          % stats of the center segment separated in half
          opts = metamerOpts(zeros(imageSize), 'windowType=flexible', 'verbose=0', 'Na=7');
          opts.flexibleMasks(1,:,:) = splitMask{h};
          masks = mkFlexibleMasks(opts);
          paramsTemp = metamerAnalysisCorr(double(imageBSD), masks, opts);
          paramsTemp.oim = NA;
          paramsNeighborHalf(nhi) = paramsTemp; 
        end
      end
    end
  end 
end

% tidy up the computed parameters
paramsCenter = match_param_struct(paramsCenter);
[paramsCenter parsNames] = pars2mat_metamers(paramsCenter);
centerSegment(1:length(centerImage)) = {0};
centerType(1:length(centerImage)) = {'center'};

paramsCenterHalf = match_param_struct(paramsCenterHalf);
[paramsCenterHalf parsNames] = pars2mat_metamers(paramsCenterHalf);
halfCenterType(1:length(halfCenterImage)) = {'half_center'};

paramsNeighbor = match_param_struct(paramsNeighbor);
[paramsNeighbor parsNames] = pars2mat_metamers(paramsNeighbor);
neighborType(1:length(neighborImage)) = {'neighbor'};

paramsNeighborHalf = match_param_struct(paramsNeighborHalf);
[paramsNeighborHalf parsNames] = pars2mat_metamers(paramsNeighborHalf);
halfNeighborType(1:length(halfNeighborImage)) = {'half_neighbor'};

% put everything together in one cell
params = [parsNames; num2cell(paramsCenter); num2cell(paramsCenterHalf); ...
  num2cell(paramsNeighbor); num2cell(paramsNeighborHalf)];
imageCol = ['ImageName'; centerImage'; halfCenterImage'; neighborImage'; halfNeighborImage'];
typeCol = ['Type'; centerType'; halfCenterType'; neighborType'; halfNeighborType'];
segmentCol = ['Segment'; centerSegment'; halfCenterNumber'; neighborSegment'; halfNeighborSegment'];
exportCell = [imageCol, typeCol, segmentCol, params];

cell2csv("BSD_stats_Corr.csv", exportCell, savingDir);

