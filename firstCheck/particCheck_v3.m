
clearvars -except sXlDat sXl ; close all;
addpath('../func');

% sim parameters
extractData             = 1;
sXl.participantID       = '327';
vFeatures = {'PupilRight', 'PupilLeft', 'GazePointX', 'GazePointY'};
sCsf.featureName        = 'GazePointY' ;  % options are: PupilRight, PupilLeft, GazePointX, GazePointY .
sCsf.eigVecNum          = 4;
%sCsf.type               = 'pca';
sCsf.type               = 'diffusion maps';
sCsf.sDifMap.gausFactor = 8;
sCsf.minVecLength       = 40;
sCsf.maxGazePointY      = 500;
sCsf.minGazePointY      = 200;

%%  set constants, extract data from excel by feature, get lines idx
if extractData==1
    [sXl,sXlDat] = xlExtract('../../Data/Ordered/', sXl.participantID);
end

%% get feature data splitted to lines
for i=1:length(vFeatures)
    sDat.(vFeatures{i}) = split2vecs(sXlDat.(vFeatures{i}),sXlDat.mLinesIdx(:,1),sXlDat.mLinesIdx(:,2));
end
%%  filter outliers by gazepointY 
sDat.GazePointYoriginal = sDat.GazePointY;
sDat.GazePointXoriginal = sDat.GazePointX;
sDat.PupilLeftoriginal  = sDat.PupilLeft;
sDat.PupilRightoriginal = sDat.PupilRight;
% sDat.GazePointX.M(sDat.GazePointY.M > sCsf.maxGazePointY | sDat.GazePointY.M < sCsf.minGazePointY) = nan;
% sDat.PupilLeft.M(sDat.GazePointY.M > sCsf.maxGazePointY | sDat.GazePointY.M < sCsf.minGazePointY) = nan;
% sDat.PupilRight.M(sDat.GazePointY.M > sCsf.maxGazePointY | sDat.GazePointY.M < sCsf.minGazePointY) = nan;
% sDat.GazePointY.M(sDat.GazePointY.M > sCsf.maxGazePointY | sDat.GazePointY.M < sCsf.minGazePointY) = nan;


%%%  start manipulating data
%% plot eye gaze point x-y

close all;
for i=1:sXlDat.numOfLines
    tmpX = sDat.GazePointX.M(:,i);
    tmpY = sDat.GazePointY.M(:,i);
%      indStop(i) = find(sDat.GazePointX.M(:,i) > min(tmpX) + 0.9*(max(tmpX)-min(tmpX)),1);
%     while ~isempty(tmpX)
%         [~,indStart(i)] = min(tmpX);
%         if indStart(i) < indStop(i)
%             break;
%         else
%             tmpX(indStart(i))= [];
%         end
%     end
%     tmpX = sDat.GazePointX.M(indStart(i):indStop(i),i);
%     tmpY = sDat.GazePointY.M(indStart(i):indStop(i),i); 
    tmpX(isnan(tmpX))=[]; tmpY(isnan(tmpY))=[];
    if mod(i,12)==1
        figure;
    end
    if mod(i,12)==0 ,subplot(3,4,12);
    else subplot(3,4,mod(i,12));
    end
    scatter(tmpX,tmpY,10,(1:length(tmpX)),'Fill');
    if mod(i,12)==1
        colorbar;title(strcat('line',{' '},num2str(i), ' gazepoint'));
    else
        title(strcat('line',{' '},num2str(i)));
    end
    xlim([250 1000]); ylim([200 600]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);
end

%% plot basic lines comparison
for j = 1:length(vFeatures)
    for i=1:sXlDat.numOfLines       
         tmp =  sDat.(vFeatures{j}).M(:,i); 
%         tmp =  sDat.(vFeatures{j}).M(indStart(i):indStop(i),i);
        tmp(isnan(tmp))  = [];
        sOut.(vFeatures{j}).std(i)    = std(tmp);
        sOut.(vFeatures{j}).mean(i)   = mean(tmp);
        % diff
        
        tmpDiff = diff(tmp);
        sOut.(vFeatures{j}).negDiffCnt(i) = sum(tmpDiff<0);
        sOut.(vFeatures{j}).diffStd(i)   = std(tmpDiff);
        sOut.(vFeatures{j}).diffMean(i)  = mean(tmpDiff);
    end
       figure; subplot(3,2,1); plot(sOut.(vFeatures{j}).mean);      title(strcat('Mean',{'   '},sXl.participantID,{'-'},vFeatures{j}));
               subplot(3,2,2); plot(sOut.(vFeatures{j}).std);       title('STD');
               subplot(3,2,3); plot(sOut.(vFeatures{j}).diffMean);  title('Diff Mean');
               subplot(3,2,4); plot(sOut.(vFeatures{j}).diffStd);   title('Diff Std');   
               subplot(3,2,5); plot(sOut.(vFeatures{j}).negDiffCnt);   title('negative diff cnt'); 
end


%%

%% fft
% if sDat.minVecLength < sCsf.minVecLength;
%     sSpec.minLength = sCsf.minVecLength;
% else
%     sSpec.minLength = sDat.minVecLength;
% end
% sSpec.M = zeros(sXlDat.numOfSamp,sSpec.minLength);
% for i=1:sXlDat.numOfSamp
%     tmp  = sDat.M(~isnan(sDat.M(:,i)),i);
%     if length(tmp) < sSpec.minLength
%         sSpec.M(i,:) = zeros(1,sSpec.minLength);
%     else
%         sSpec.M(i,:) = abs(fft(tmp,sSpec.minLength));
%     end
% end

%  %%
% sCsf.M = sSpec.M;
% if strcmp(sCsf.type,'diffusion maps')
%     [mPhi, mLam] = gausDiffusinMaps(sCsf.M,'median',sCsf.sDifMap.gausFactor,sCsf.eigVecNum);
% elseif strcmp(sCsf.type,'pca')
%     [mPhi, mLam] = pcaWrap(sCsf.M,1,sCsf.eigVecNum);
% end

%%
%close all;
% %% spec plot
% sPlot.idStr = [sXl.participantID '   fft   ' sCsf.featureName];
% sPlot.pointSizfigure; plot(sSpec.M(5,:)); hold on; plot(sSpec.M(6,:)); title({sPlot.idStr '2 Sentences FFTs'});
% e = 30; sPlot.fontSize = 8;
% figure; scatter3(mPhi(:,1), mPhi(:,2), mPhi(:,3), sPlot.pointSize, 1:size(sCsf.M,1), 'Fill'); colorbar;
% sPlot.lineNums = cellfun(@num2str, num2cell(1:size(sCsf.M,1)), 'UniformOutput', false);
% text(mPhi(:,1), mPhi(:,2), mPhi(:,3),sPlot.lineNums,'FontSize',sPlot.fontSize,'Color','red'); 
% title({sPlot.idStr 'First 3 Eigen Vectors'} );
% figure; scatter3(mPhi(:,2), mPhi(:,3), mPhi(:,4), sPlot.pointSize, 1:size(sCsf.M,1), 'Fill'); colorbar;
% text(mPhi(:,2), mPhi(:,3), mPhi(:,4),sPlot.lineNums,'FontSize',sPlot.fontSize,'Color','red'); 
% title({sPlot.idStr '2,3,4 Eigen Vectors'} );
% % figure; plot(mPhi(:,1));
% % figure; plot(mPhi(:,2));
% % figure; plot(mPhi(:,3));

clear('tmp','i');