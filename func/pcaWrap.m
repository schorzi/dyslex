function [ mPhi,scores,eigenValues ] = pcaWrap( M, normalizeData,princCompNum )
if normalizeData == 1
    M = bsxfun(@rdivide,M,std(M));
end
if exist('princCompNum','var')
[mPhi,scores,eigenValues] = pca(M,'NumComponents',princCompNum);
else
[mPhi,scores,eigenValues] = pca(M);
end
mPhi = M*mPhi;

end

