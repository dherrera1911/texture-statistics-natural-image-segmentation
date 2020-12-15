##### Generate statistics data from texture database ######
pkg load image
addpath(genpath("./functions"));

texturesDir = "../../images/experiment_textures/";
savingDir = "../../data/experiment_stats/";

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

textureNames = read_unhidden(texturesDir);

paramsOrig = struct();
paramsMatched = struct();
namesCol = {};

for imN = 1:length(textureNames)
  tName = textureNames{imN};
  original = imread([texturesDir, tName, "/original.tif"]);
  matched = imread([texturesDir, tName, "/matched.tif"]);
  paramsOrig(imN) = textureAnalysisCorr(double(original), Nsc, Nor, Na);
  paramsMatched(imN) = textureAnalysisCorr(double(matched), Nsc, Nor, Na);
  namesCol{imN} = tName;
end


[parsMatOrig parsNamesOrig] = pars2mat_texture(paramsOrig);
[parsMatMatched parsNamesMatched] = pars2mat_texture(paramsMatched);
statsMat = num2cell([parsMatOrig; parsMatMatched]);
statsCell = [parsNamesOrig; statsMat];

kindCol = [repmat({"original"}, 4, 1); repmat({"matched"}, 4, 1)];
kindCell = {"type", kindCol{:}}';
textureCell = {"texture", namesCol{:}, namesCol{:}}';

exportCell = [textureCell, kindCell, statsCell];

cell2csv("experiment_textures_stats.csv", exportCell, savingDir);

