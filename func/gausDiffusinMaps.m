function [ mPhi,mLam ] = gausDiffusinMaps(M,epsType,epsFactor,eigenVecNum )

mW  = squareform( pdist(M) );
if strcmp(epsType,'median')
  eps = epsFactor * median(mW(:)); 
end 
mK  = exp( -mW.^2 / eps.^2 );
mA  = bsxfun(@rdivide, mK, sum(mK, 2));
% mD  = diag(sum(mK));
% mA  = mD^(-1/2) * mK * mD^(-1/2);
if exist('eigenVecNum','var');
    [mPhi, mLam] = eigs(mA, eigenVecNum+1);
else
    [mPhi, mLam] = eigs(mA);
end
mPhi  = mPhi(:,2:end);
mLam  = diag(mLam(2:end,2:end));

end

