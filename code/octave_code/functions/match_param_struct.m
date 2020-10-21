%{
Transform structure of the params struct obtained from calculating for
several different images (as done within the get_stimulus_params
script) into the structure returned by the get_target_params script
(mainly removes the 'mask' layer)

inStruct is a struct array, with only one mask used for calculating
the params
%}

function matchedStruct = match_param_struct(inStr)
  function out = remove_mask(in)
    out = cellfun(@(x) x.mask{1}, in, 'UniformOutput', false);
  end
  nScales = length(inStr(1).autoCorrMag.scale);
  for t = 1:length(inStr)
    inStr(t).LPskew = [inStr(t).LPskew.scale{:}];
    inStr(t).LPkurt = [inStr(t).LPkurt.scale{:}];
    inStr(t).autoCorrReal.scale = remove_mask(inStr(t).autoCorrReal.scale);
    inStr(t).autoCorrReal.vHPR = inStr(t).autoCorrReal.vHPR.mask{1};
    inStr(t).autoCorrRealFull.scale = remove_mask(inStr(t).autoCorrRealFull.scale);
    for s = 1:nScales
      inStr(t).autoCorrMag.scale{s}.ori = ...
        remove_mask(inStr(t).autoCorrMag.scale{s}.ori);
      inStr(t).autoCorrMagFull.scale{s}.ori = ...
        remove_mask(inStr(t).autoCorrMagFull.scale{s}.ori);
    end
    inStr(t).magMeans = [inStr(t).magMeans.band{:}]'; 
    inStr(t).cousinMagCorr.scale = remove_mask(inStr(t).cousinMagCorr.scale);
    inStr(t).parentMagCorr.scale = remove_mask(inStr(t).parentMagCorr.scale);
    inStr(t).cousinRealCorr.scale = remove_mask(inStr(t).cousinRealCorr.scale);
    inStr(t).parentRealCorr.scale = remove_mask(inStr(t).parentRealCorr.scale);
  end
  matchedStruct = inStr;
end


