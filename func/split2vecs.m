function [ mat_out,minVecLength,maxVecLength ] = split2vecs(vec,startInd,endInd)

deltas         = (endInd - startInd)+1;
maxVecLength   = max(deltas);
minVecLength = min(deltas);
tmp            = zeros(maxVecLength,length(startInd));

for i=1:length(startInd)
        tmp(1:deltas(i),i)     = vec(startInd(i):endInd(i));
        tmp(deltas(i)+1:end,i) = nan;
end

mat_out = tmp;

end

