function [ out ] = partOfVector(vec,startInd,endInd)
tmp = [];
for i=1:length(startInd)
        tmp = [tmp;vec(startInd(i):endInd(i))];
end

out = tmp;

end

