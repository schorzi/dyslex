function [ sDat ] = split2vecs(vec,startInd,endInd)

maxLength     = max(endInd - startInd);
tmp           = zeros(maxLength,length(startInd));
vecLength     = zeros(length(startInd),1);
for i=1:length(startInd)
        vTmp                   = vec(startInd(i):endInd(i));
        vTmp(isnan(vTmp))      = [];
        vecLength(i)           = length(vTmp);
        tmp(1:vecLength(i),i)  = vTmp;
        tmp(vecLength(i)+1:end,i) = nan;
end

sDat.maxVecLength  = max(vecLength);
sDat.minVecLength  = min(vecLength);
sDat.meanVecLength = mean(vecLength);
sDat.medVecLenght  = median(vecLength);
sDat.M = tmp;

end

