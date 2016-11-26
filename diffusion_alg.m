close all;
clear;

D = 10;
N = 500;

mX1 = randn(D, N) + 1;
mX2 = randn(D, N) - 1;
mX  = [mX1, mX2];
vY  = [zeros(N, 1); ones(N, 1)];

mW  = squareform( pdist(mX') );
eps = .5 * median(mW(:).^2);
mK  = exp(- mW.^2 / eps);
mA  = bsxfun(@rdivide, mK, sum(mK, 2));

[mU, mL, mV] = eig(mA);

mP1 = mU(:, 2);
mP2 = mU(:, 2:3);
mP3 = mU(:, 2:4);

figure; plot(mP1);
figure; scatter(mP2(:,1), mP2(:,2), 100, mW(1,:)', 'Fill');       colorbar;
figure; scatter3(mP3(:,1), mP3(:,2), mP3(:,3), 100, vY', 'Fill'); colorbar;


figure; imagesc(mW); colorbar;
figure; imagesc(mK); colorbar;

