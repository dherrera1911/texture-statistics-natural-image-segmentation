%{
Turn a struct array containing the statistical parameters
of the textures into a matrix, and a cell array containing
the names of the parameters

Removes some parameters as indicated by:
scalesOut = vector with numbers indicating the scales to 
leave out of the matrix
nbhSize = size of the neighborhood to keep for the correlations
through space

This function doesn't include the parameters for the highpass

%}

function [parsMat parsNames] = pars2mat_metamers(params, scalesOut, nbhSize)
  nTrials = length(params); 
  nScales = size(params(1).autoCorrMag.scale, 2);
  nOri = size(params(1).autoCorrMag.scale{1}.ori,2);

  magMeansOut = [];
  switch nargin
    case 1
      scalesOut = [];
      nbhPruning = false;
      nbhSize = size(params(1).autoCorrReal.scale{1},1);
      remainingScales = 1:nScales;
    case 2
      nbhPruning = false;
      nbhSize = size(params(1).autoCorrReal.scale{1},1);
      for j = 1:length(scalesOut)
        bandI = nOri*(scalesOut(j)-1)+1;
        bandF = nOri*scalesOut(j);
        magMeansOut = [magMeansOut bandI:bandF];
      end
      remainingScales = setdiff(1:nScales, scalesOut);
    case 3
      nbhPruning = true;
      for j = 1:length(scalesOut)
        bandI = nOri*(scalesOut(j)-1)+1;
        bandF = nOri*scalesOut(j);
        magMeansOut = [magMeansOut bandI:bandF];
      end
      origNbh = size(params(1).autoCorrReal.scale{1},1);
      nbhI = [origNbh-nbhSize]/2+1;
      nbhF = nbhI+nbhSize-1;
      remainingScales = setdiff(1:nScales, scalesOut);
  end

  nanCorr = tril(ones(nbhSize));  
  nanCorr(nanCorr==0) = NaN;
  nanAutoCorr = nanCorr;
  for i = nbhSize:-1:[(nbhSize+1)/2+1]
    nanAutoCorr(i,i) = NaN;
  end
  nanOriCorr = 1-triu(ones(nOri));
  nanOriCorr(nanOriCorr==0) = NaN;

  pixStatsMat = [params.pixelStats]';
  LPskewMat = [params.LPskew];
  LPskewMat = reshape(LPskewMat, nScales+1, nTrials);
  LPskewMat = LPskewMat';
  LPskewMat(:,scalesOut) = [];
  LPkurtMat = [params.LPkurt];
  LPkurtMat = reshape(LPkurtMat, nScales+1, nTrials);
  LPkurtMat = LPkurtMat';
  LPkurtMat(:,scalesOut) = [];

  acrMat = [];
  acmMat = [];
  mmMat = [];
  cmcMat = [];
  pmcMat = [];
  crcMat = [];
  prcMat = [];

  for t = 1:length(params)
    %remove unnecessary matrices and zeros
    params(t).cousinRealCorr.scale(nScales+1) = [];%uninterpretable
    params(t).parentRealCorr.scale(nScales+1) = [];%empty
    params(t).parentRealCorr.scale(nScales) = [];%uninterpretable
    params(t).parentMagCorr.scale(nScales:nScales+1) = [];%empty
    params(t).cousinRealCorr.scale(1:end) = cellfun(@(x) x(1:4, 1:4), ...
      params(t).cousinRealCorr.scale(1:end), 'UniformOutput', false);
    params(t).parentRealCorr.scale(1:end) = cellfun(@(x) x(1:4, 1:8), ...
      params(t).parentRealCorr.scale(1:end), 'UniformOutput', false);

    %remove scales solicited for removal
    params(t).autoCorrReal.scale(scalesOut) = [];
    params(t).autoCorrMag.scale(scalesOut) = [];
    params(t).magMeans = params(t).magMeans(2:end);
    params(t).magMeans(magMeansOut) = [];
    params(t).cousinMagCorr.scale(scalesOut) = [];
    params(t).parentMagCorr.scale(scalesOut) = [];
    params(t).cousinRealCorr.scale(scalesOut) = [];
    params(t).parentRealCorr.scale(scalesOut) = [];

    %if solicited, prune the correlation neighborhoods
    if nbhPruning
      params(t).autoCorrReal.scale = cellfun(@(x) x(nbhI:nbhF,nbhI:nbhF), ...
        params(t).autoCorrReal.scale, 'UniformOutput', false);
      for s = 1:length(params(t).autoCorrMag.scale)
        params(t).autoCorrMag.scale{s}.ori = cellfun(@(x) x(nbhI:nbhF,nbhI:nbhF), ...
          params(t).autoCorrMag.scale{s}.ori, 'UniformOutput', false);
      end
    end
     
    %change values repeated in correlations for NaNs
    params(t).autoCorrReal.scale = cellfun(@(x) x.*nanAutoCorr, ...
      params(t).autoCorrReal.scale, 'UniformOutput', false);
    for s = 1:length(params(t).autoCorrMag.scale)
      params(t).autoCorrMag.scale{s}.ori = cellfun(@(x) x.*nanAutoCorr, ...
        params(t).autoCorrMag.scale{s}.ori, 'UniformOutput', false);
    end
    params(t).cousinMagCorr.scale = cellfun(@(x) x.*nanOriCorr, ...
      params(t).cousinMagCorr.scale, 'UniformOutput', false);
    params(t).cousinRealCorr.scale = cellfun(@(x) x.*nanOriCorr, ...
      params(t).cousinRealCorr.scale, 'UniformOutput', false);
     
    %convert the remaining values into matrices
    acrRow = [params(t).autoCorrReal.scale{:}];
    acrMat(t,:) = acrRow(:);
    acmRow = [];
    for s = 1:length(params(t).autoCorrMag.scale)
      acmRow = [acmRow params(t).autoCorrMag.scale{s}.ori{:}];
    end
    acmMat(t,:) = acmRow(:);
    mmMat(t,:) = params(t).magMeans;
    cmcRow = [params(t).cousinMagCorr.scale{:}];
    cmcMat(t,:) = cmcRow(:); 
    pmcRow = [params(t).parentMagCorr.scale{:}];
    pmcMat(t,:) = pmcRow(:);
    crcRow = [params(t).cousinRealCorr.scale{:}];
    crcMat(t,:) = crcRow(:);
    prcRow = [params(t).parentRealCorr.scale{:}];
    prcMat(t,:) = prcRow(:);
  end
  

  %make names of columns
  scalesCell = num2cell(remainingScales);
  scalesCellLP = {scalesCell{:}, 'LP'};
  oriCell = num2cell(1:nOri);
  nbhNums = num2cell(-(nbhSize-1)/2:(nbhSize-1)/2);
  parentsStatsScalesNames = {"S1S2", "S2S3", "S3S4", "S4S5", "S5S6"}; 
  parentsStatsScalesNames = parentsStatsScalesNames(remainingScales(1:end-1));

  statsName = {'pix_mean', 'pix_var', 'pix_skew', 'pix_kurt', 'pix_min', 'pix_max'};
  skewName = all_combs({'LPskew_S'}, scalesCellLP);
  kurtName = all_combs({'LPkurt_S'}, scalesCellLP);
  acrName = all_combs({'acr_S'}, scalesCellLP, {'X'}, nbhNums, {'Y'}, nbhNums); 
  acmName = all_combs({'acm_S'}, scalesCell, {'O'}, oriCell, {'X'}, ...
    nbhNums, {'Y'}, nbhNums);
  mmName = all_combs({'mm_S'}, scalesCell, {'O'}, oriCell);
  mmName{end+1} = 'mm_LP';
  cmcName = all_combs({'cmc_S'}, scalesCell, {'O'}, oriCell, {'O'}, oriCell);
  pmcName = all_combs({'pmc_'}, parentsStatsScalesNames, {'O'}, oriCell, {'O'}, oriCell);
#  crcName = all_combs({'crc_S'}, scalesCell, {'O'}, oriCell, {'O'}, oriCell);
  prcName = all_combs({'prc_'}, parentsStatsScalesNames, {'Ph'}, num2cell(1:2), {'O'}, ...
    oriCell, {'O'}, oriCell);

  %Append the different param matrices and the names cells in the
  %same order. Eliminate the NaN columns
  parsMat = [pixStatsMat, LPskewMat, LPkurtMat, acrMat, acmMat, ...
    mmMat, cmcMat, pmcMat, prcMat];
  parsNames = {statsName{:}, skewName{:}, kurtName{:}, acrName{:}, ...
    acmName{:}, mmName{:}, cmcName{:}, pmcName{:}, prcName{:}};

  removedValues = isnan(parsMat(1,:));  
  parsMat = parsMat(:,~removedValues);
  parsNames = parsNames(~removedValues);

end

