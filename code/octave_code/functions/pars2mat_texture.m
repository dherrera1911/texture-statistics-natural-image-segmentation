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

function [parsMat parsNames] = pars2mat_texture(params)

  nSamples = length(params); 
  nScales = size(params(1).autoCorrMag, 3);
  nOri = size(params(1).autoCorrMag, 4);
  nbhSize = size(params(1).autoCorrMag, 1);

  nanCorr = tril(ones(nbhSize));  
  nanCorr(nanCorr==0) = NaN;
  nanCousinCorr = tril(ones(nOri));  
  nanCousinCorr(nanCousinCorr==0) = NaN;
  nanAutoCorr = nanCorr;
  # redundant indicator for orientations
  nanOriCorr = 1-triu(ones(nOri));
  nanOriCorr(nanOriCorr==0) = NaN;

  # remove redundant numbers in the diagonal
  for i = nbhSize:-1:[(nbhSize+1)/2+1]
    nanAutoCorr(i,i) = NaN;
  end
  # replicate matrix for across scales, orientations
  nanAutoCorrReal = repmat(nanAutoCorr, 1, 1, nScales+1);
  nanAutoCorrMag = repmat(nanAutoCorr, 1, 1, nScales, nOri);
  nanOriCorr = repmat(nanOriCorr, 1, 1, nScales);


  pixelStatsMat = [];
  LPskewMat = [];
  LPkurtMat = [];
  acrMat = [];
  acmMat = [];
  mmMat = [];
  cmcMat = [];
  pmcMat = [];
  crcMat = [];
  prcMat = [];

  for t = 1:length(params)

    pixelStatsMat = [pixelStatsMat; params(t).pixelStats];
    LPskewMat = [LPskewMat; params(t).pixelLPStats(:,1)'];
    LPkurtMat = [LPkurtMat; params(t).pixelLPStats(:,2)'];

    %remove unnecessary matrices and zeros
    params(t).magMeans = params(t).magMeans(2:end);
    params(t).cousinRealCorr(:,:,nScales+1) = [];%uninterpretable
    params(t).parentRealCorr(:,:,nScales) = [];%uninterpretable
    params(t).parentMagCorr(:,:,nScales) = [];%empty
    params(t).cousinMagCorr(:,:,nScales+1) = [];%empty

    params(t).cousinRealCorr = params(t).cousinRealCorr(1:nOri,1:nOri,:);
    params(t).parentRealCorr = params(t).parentRealCorr(1:nOri,1:(2*nOri),:);

    %change values repeated in correlations for NaNs
    params(t).autoCorrReal = params(t).autoCorrReal .* nanAutoCorrReal;
    params(t).autoCorrMag = params(t).autoCorrMag .* nanAutoCorrMag;
    params(t).cousinMagCorr = params(t).cousinMagCorr .* nanOriCorr;
    params(t).cousinRealCorr = params(t).cousinRealCorr .* nanOriCorr;

    %convert the remaining values into matrices
    acrMat(t,:) = params(t).autoCorrReal(:);
    acmMat(t,:) = params(t).autoCorrMag(:);
    mmMat(t,:) = params(t).magMeans;
    cmcMat(t,:) = params(t).cousinMagCorr(:);
    pmcMat(t,:) = params(t).parentMagCorr(:);
    crcMat(t,:) = params(t).cousinRealCorr(:);
    prcMat(t,:) = params(t).parentRealCorr(:);
  end

  %make names of columns
  scalesCell = num2cell(1:nScales);
  scalesCellLP = {scalesCell{:}, 'LP'};
  oriCell = num2cell(1:nOri);
  nbhNums = num2cell(-(nbhSize-1)/2:(nbhSize-1)/2);
  parentsStatsScalesNames = {"S1S2", "S2S3", "S3S4", "S4S5", "S5S6"}; 
  parentsStatsScalesNames = parentsStatsScalesNames(1:(nScales-1));

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
  parsMat = [pixelStatsMat, LPskewMat, LPkurtMat, acrMat, acmMat, ...
    mmMat, cmcMat, pmcMat, prcMat];
  parsNames = {statsName{:}, skewName{:}, kurtName{:}, acrName{:}, ...
    acmName{:}, mmName{:}, cmcName{:}, pmcName{:}, prcName{:}};

  removedValues = isnan(parsMat(1,:));  
  parsMat = parsMat(:,~removedValues);
  parsNames = parsNames(~removedValues);

end

