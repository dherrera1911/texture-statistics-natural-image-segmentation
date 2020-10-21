%%%%%%% Generate statistics data from texture database %%%%%%%%%
pkg load image

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

%{
Get a segmentation map as input, and initial coords for a given pixel.
Return a vector indicating the neighboring segments up, down, left and
right, returning NA if no neighboring segment in that direction.
initialCoords are in [row col]
%}
function neighborSegments = find_neighbor(segmentMap, initialCoords)
  xi = initialCoords(2);
  yi = initialCoords(1);
  mainSeg = segmentMap(yi, xi);
  % find neighbouring segments
  left = first_different(flip(segmentMap(yi, 1:xi)), mainSeg);
  right = first_different(segmentMap(yi, xi:end), mainSeg);
  up = first_different(flip(segmentMap(1:yi, xi)), mainSeg);
  down = first_different(segmentMap(yi:end, xi), mainSeg);
  neighborSegments = [left, right, up, down];
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

% return first value of vector different from reference value
function outVal = first_different(inputVec, referenceVal)
  ind = find(inputVec != referenceVal);
  if isempty(ind)
    outVal = NA;
  else
    outVal = double(inputVec(ind(1)));
  end
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

%{
Normalize acrs
%}
function normAcr = normalize_acr(inputAcr)

end
  

%%%%%%%%%%%% Compute the stats %%%%%%%%%%%%%%%%
imageFiles = read_unhidden(imageDir);
segmentFiles = read_unhidden(segmentsDir);

minSize = 64*64;

im = 2;
imName = strsplit(imageFiles{im}, "."){1};
% load image and segmentation map
imageBSD = imread([imageDir, imageFiles{im}]);
imageBSD = rgb2gray(imageBSD);

segment = load([segmentsDir, segmentFiles{im}]);
segmentMap = segment.groundTruth{2}.Segmentation;

%get sizes and center coords
imageSize = size(imageBSD);
xCenter = round(imageSize(2)/2);
yCenter = round(imageSize(1)/2);

% get center segment and neighbor segments
centerSeg = segmentMap(yCenter, xCenter);
neighborSegments = unique(find_neighbor(segmentMap, [yCenter, xCenter]));
neighborSegments = neighborSegments(!isna(neighborSegments)); 
splitMask = split_mask(centerMask);

% filter by minimum size
neighborSegments = filter_segments(neighborSegments, segmentMap, minSize);
centerMask = double(segmentMap == centerSeg);
sizeCenter = sum(centerMask(:));

imwrite(imageBSD, "fullImage.png")
imwrite(mat2gray(segmentMap), "segmentMap.png")
imwrite(imageBSD, "centerFull.png", 'Alpha', double(segmentMap==34)*255)
imwrite(imageBSD, "centerLeft.png", 'Alpha', splitMask{1}*255)
imwrite(imageBSD, "centerRight.png", 'Alpha', splitMask{2}*255)
imwrite(imageBSD, "neighbor.png", 'Alpha', double(segmentMap==8)*255)




### export texture patches ###
texturesDir = "../../images/textures_initial/";
addpath(genpath("../octave_code/functions"));

# read contents of a dir that are not hidden files
function [files] = read_unhidden(inputDir)
  temp = dir(inputDir);
  temp = {temp.name};
  hiddenFiles = regexp(temp, "^\\.");
  unhiddenFiles = cellfun("isempty", hiddenFiles);
  files = temp(unhiddenFiles);
end

# return whether the strings in a cell are an image file type
function [isImage] = is_image_file(strCell)
  imageTypes = {".tif", ".jpg", ".png", ".tiff", ".pnm"};
  isImage = zeros(length(strCell), 1);
  for it = 1:length(imageTypes)
    matchInd = strfind(lala, imageTypes{it});
    matchInd = !(cellfun("isempty", matchInd));
    isImage = isImage + matchInd;
  end
end

function [xi, yi, xf, yf] = quadrant_coords(quadrant, imSize, quadrantSize)
  switch quadrant
    case 1
      xf = imSize(2);
      xi = xf - quadrantSize + 1; 
      yi = 1;
      yf = yi + quadrantSize - 1; 
    case 2
      xi = 1; 
      xf = xi + quadrantSize - 1;
      yi = 1;
      yf = yi + quadrantSize - 1; 
    case 3
      xi = 1; 
      xf = xi + quadrantSize - 1;
      yf = imSize(1); 
      yi = yf - quadrantSize + 1;
    case 4
      xf = imSize(2);
      xi = xf - quadrantSize + 1; 
      yf = imSize(1); 
      yi = yf - quadrantSize + 1;
  end
end

textureFiles = read_unhidden(texturesDir);
quadrantSize = 128;

t = 0;
imN = 7;

texture = imread([texturesDir, textureFiles{imN}]);
if (ndims(texture) == 3)
  texture = rgb2gray(texture);
end
grayTexture = mat2gray(texture);

% set mean-var of texture
normalizedStatsTexture = (grayTexture - mean(grayTexture(:))) / std(grayTexture(:));
normalizedStatsTexture = uint8(((normalizedStatsTexture * 0.2) + 0.5)*255);

%textureName = strsplit(textureFiles{imN}, "."){1};
for quadrant = 1:4
  t = t + 1;
  [xi, yi, xf, yf] = quadrant_coords(quadrant, size(texture), quadrantSize);
  quadrantCol{t} = quadrant;
  % sample and compute statistics
  textureQuadrantStats{quadrant} = normalizedStatsTexture(yi:yf, xi:xf);
end
imshow(normalizedStatsTexture)

imwrite(textureQuadrantStats{1}, "textureSameQ1.png")
imwrite(textureQuadrantStats{2}, "textureSameQ2.png")

imwrite(textureQuadrantStats{1}, "textureDiffQ1.png")
imwrite(textureQuadrantStats{2}, "textureDiffQ2.png")




