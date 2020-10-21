##### Generate statistics data from texture database ######
pkg load image
addpath(genpath("./functions"));

texturesDir = "../../images/textures_initial/";
savingDir = "../../data/texture_stats/";

Nsc = 4;
Nor = 4;
Na = 7;

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

# build a tree with the data for the textures
#function [fileTreeStruct, imagePaths] = recursive_image_search(inputDir)
#  fileTreeStruct = struct("dirName", inputDir, "children", []);
#  dirContents = read_unhidden(inputDir);
#  
#  isDir = cellfun("isdir", strcat(inputDir, dirContents));
#  isImage = is_image_file(dirContents);
#end
   
textureFiles = read_unhidden(texturesDir);
quadrantSize = 128;
params = struct();
paramsPixNorm = struct();
paramsStatsNorm = struct();
namesCol = {};
quadrantCol = {};

t = 0;
for imN = 1:length(textureFiles)
  texture = imread([texturesDir, textureFiles{imN}]);
  if (ndims(texture) == 3)
    texture = rgb2gray(texture);
  end
  grayTexture = mat2gray(texture);
  % set maximum min of texture
  normalizedPixTexture = 255*(grayTexture);
  % set mean-var of texture
  normalizedStatsTexture = (grayTexture - mean(grayTexture(:))) / std(grayTexture(:));
  normalizedStatsTexture = uint8(((normalizedStatsTexture * 0.2) + 0.5)*255);
  %textureName = strsplit(textureFiles{imN}, "."){1};
  textureName = textureFiles{imN}(1:(end-4));
  for quadrant = 1:4
    t = t + 1;
    [xi, yi, xf, yf] = quadrant_coords(quadrant, size(texture), quadrantSize);
    namesCol{t} = textureName;
    quadrantCol{t} = quadrant;
    % sample and compute statistics
    textureQuadrant = texture(yi:yf, xi:xf);
    %params(t) = textureAnalysisCorr(double(textureQuadrant), Nsc, Nor, Na);
    % compute statistics for version with normalized max-min
    %textureQuadrantPix = normalizedPixTexture(yi:yf, xi:xf);
    %paramsPixNorm(t) = textureAnalysisCorr(double(textureQuadrantPix), Nsc, Nor, Na);
    % compute statistics for version with normalized mean-var
    textureQuadrantStats = normalizedStatsTexture(yi:yf, xi:xf);
    paramsStatsNorm(t) = textureAnalysisCorr(double(textureQuadrantStats), Nsc, Nor, Na);
  end
end

quadrantCell = {"quadrant", quadrantCol{:}}';
textureCell = {"texture", namesCol{:}}';
  
# tidy SS extracted
%[parsMat parsNames] = pars2mat_texture(params);
%[parsMatPixNorm parsNamesPixNorm] = pars2mat_texture(paramsPixNorm);
[parsMatStatsNorm parsNamesStatsNorm] = pars2mat_texture(paramsStatsNorm);

%parsCell = [parsNames; num2cell(parsMat)];
%parsCellPixNorm = [parsNamesPixNorm; num2cell(parsMatPixNorm)];
parsCellStatsNorm = [parsNamesStatsNorm; num2cell(parsMatStatsNorm)];

%exportCell = [textureCell, quadrantCell, parsCell];
%exportCellPixNorm = [textureCell, quadrantCell, parsCellPixNorm];
exportCellStatsNorm = [textureCell, quadrantCell, parsCellStatsNorm];

%cell2csv("texture_stats.csv", exportCell, savingDir);
%cell2csv("texture_stats_pixNorm.csv", exportCellPixNorm, savingDir);
cell2csv("texture_stats_statsNorm.csv", exportCellStatsNorm, savingDir);

