close all;
clear;

path     = 'C:\Users\Oryair\Desktop\Workarea\recorded\1\';
fileName = 'sentences_new_Data_Export_pair with audio.xlsx';
% data     = xlsread([path, fileName], 'BI:BT');
data     = xlsread([path, fileName], 'BY:CB');

%%
M          = data;
mNaN       = isnan(M);
vNaN       = logical(sum(mNaN, 2));
M(vNaN, :) = [];
M          = M(1:1:end,:);
%%
figure; imagesc(M); colorbar;
figure; wiggle(M); colorbar;

%%
mW  = squareform( pdist(M) );
eps = 50 * median(mW(:));
mK  = exp( -mW.^2 / eps.^2 );
mA  = bsxfun(@rdivide, mK, sum(mK, 2));
% mD  = diag(sum(mK));
% mA  = mD^(-1/2) * mK * mD^(-1/2);
[mPhi, mLam] = eigs(mA, 4);

%%
figure; scatter3(mPhi(:,2), mPhi(:,3), mPhi(:,4), 100, 1:size(M,1), 'Fill'); colorbar;