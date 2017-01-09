
clearvars -except sRawXlDat ; close all;
addpath('../func');

%% sim parameters
extractData             = 1;
sXl.participantID       = '295';
sXl.sixLineStart        = 37000;     %in micro seconds %252=19700 295=17221
sCsf.featureName        = 'GazePointY' ;  % options are: PupilRight, PupilLeft, GazePointX, GazePointY .
sCsf.eigVecNum          = 4;
%sCsf.type               = 'pca';
sCsf.type               = 'diffusion maps';
sCsf.sDifMap.gausFactor = 8;
sCsf.minVecLength       = 40;
sC.Fs                   = 44100;
%% add participant data folder
sXl.dataFolderPath           = '../../Data/Ordered/';
sXl.participantFolderPath    = strcat(sXl.dataFolderPath, sXl.participantID);
sXl.participantFilePathNoExt = strcat(sXl.participantFolderPath,'/',sXl.participantID);
sXl.xlFileName               = strcat(sXl.participantFilePathNoExt,'.xlsx');
sXl.linesIdxFile              = strcat('linesIdx',sXl.participantID);
addpath(sXl.participantFolderPath);
if extractData == 1
[sRawXlDat.num,sRawXlDat.text,sRawXlDat.raw]  = xlsread(sXl.xlFileName);
end;
%% relevent excel columns
sXl.sCol.('PupilLeft')  = 72; sXl.sCol.('PupilRight') = 73;
sXl.sCol.('GazePointX') = 57; sXl.sCol.('GazePointY') = 58; 
sXl.sCol.ValidLeft  = 74; sXl.sCol.ValidRight = 75;
sXl.sCol.EventType  = 37; sXl.sCol.EventName  = 38;
sXl.sCol.time       = 18; % in num mtx
sXl.LineStart       = 'InstructionStart'; 
sXl.LineEnd         = 'InstructionEnd';
%% get idxs of lines - only valid data
load(sXl.linesIdxFile);  % variable linesIdx
sXlDat.vTime  = sRawXlDat.num(:,sXl.sCol.time);
sXlDat.vTime(isnan(sXlDat.vTime))=[];  % time in milliseconds
sXlDat.vTime  = (sXlDat.vTime - sXlDat.vTime(1));      % align to zero
sXlDat.mLinesTime = linesIdx * (1000/sC.Fs);                % convert idxs to time in milliseconds
sXl.xl2mp3delay   = sXl.sixLineStart - sXlDat.mLinesTime(6,1);    % set mp3 to excel delay
sXlDat.mLinesTimeSync = sXlDat.mLinesTime + sXl.xl2mp3delay;
sXlDat.edges = [-Inf, mean([sXlDat.vTime(2:end)'; sXlDat.vTime(1:end-1)']), +Inf];
sXlDat.mLinesIdx = discretize(sXlDat.mLinesTimeSync, sXlDat.edges);

sXlDat.mNum   = sRawXlDat.num;
sXlDat.mTxt   = sRawXlDat.text(2:end,:);
% sXlDat.vValidEyesIdx = find((sXlDat.mNum(:,sXl.sCol.ValidRight) == 0 & sXlDat.mNum(:,sXl.sCol.ValidLeft) == 0) | ...
%                         strcmp(sXlDat.mTxt(:,sXl.sCol.EventType),sXl.LineStart) | ...
%                         strcmp(sXlDat.mTxt(:,sXl.sCol.EventType),sXl.LineEnd) );
% sXlDat.mNumValid  = sXlDat.mNum(sXlDat.vValidEyesIdx,:);
% sXlDat.mTxtValid  = sXlDat.mTxt(sXlDat.vValidEyesIdx,:);
% sXlDat.vStartInd  = find(strcmp(sXlDat.mTxtValid(:,sXl.sCol.EventType),sXl.LineStart))+1;
% sXlDat.vStartInd  = sXlDat.vStartInd(2:end);
% sXlDat.vEndInd    = find(strcmp(sXlDat.mTxtValid(:,sXl.sCol.EventType),sXl.LineEnd))-1;
% sXlDat.vEndInd    = sXlDat.vEndInd(2:end);
% if length(sXlDat.vEndInd) ~= length(sXlDat.vStartInd)
%    error('StartIdx and EndIdx must be the same size, check Xl file'); end

sXlDat.numOfSamp      = size(sXlDat.mLinesIdx,1);

%% get feature data
sXlDat.vFeatureIdx    = sXl.sCol.(sCsf.featureName);
sXlDat.vFeatureDat    = sXlDat.mNum(:,sXlDat.vFeatureIdx);
%sXlDat.vFeatureDat(isnan(sXlDat.vFeatureDat)) = [];  % remove empty
sDat = split2vecs(sXlDat.vFeatureDat,sXlDat.mLinesIdx(:,1),sXlDat.mLinesIdx(:,2));
%% spectrum of data
if sDat.minVecLength < sCsf.minVecLength;
    sSpec.minLength = sCsf.minVecLength;
else
    sSpec.minLength = sDat.minVecLength;
end
% extract
for i=1:sXlDat.numOfSamp
    tmp  = sDat.M(~isnan(sDat.M(:,i)),i);
    sOut.std(i)  = std(tmp);
    sOut.mean(i) = mean(tmp);
    tmpDiff = diff(tmp);
    sOut.vDiffMean(i) = mean(tmpDiff);
    sOut.vDiffStd(i)  = std(tmpDiff);
    sOut.vDiffBack(i) = length(find(tmpDiff<0));
end

figure; plot(sOut.std); title('STD');
figure; plot(sOut.mean); title('Mean');
figure; plot(sOut.vDiffMean); title('DiffMean');
figure; plot(sOut.vDiffStd); title('DiffStd');
% fft
sSpec.M = zeros(sXlDat.numOfSamp,sSpec.minLength);
for i=1:sXlDat.numOfSamp
    tmp  = sDat.M(~isnan(sDat.M(:,i)),i);
    if length(tmp) < sSpec.minLength
        sSpec.M(i,:) = zeros(1,sSpec.minLength);
    else
        sSpec.M(i,:) = abs(fft(tmp,sSpec.minLength));
    end
end

 %%
sCsf.M = sSpec.M;
if strcmp(sCsf.type,'diffusion maps')
    [mPhi, mLam] = gausDiffusinMaps(sCsf.M,'median',sCsf.sDifMap.gausFactor,sCsf.eigVecNum);
elseif strcmp(sCsf.type,'pca')
    [mPhi, mLam] = pcaWrap(sCsf.M,1,sCsf.eigVecNum);
end

%%
%close all;
%% spec plot
sPlot.idStr = [sXl.participantID '   fft   ' sCsf.featureName];
sPlot.pointSize = 30; sPlot.fontSize = 8;
figure; plot(sSpec.M(5,:)); hold on; plot(sSpec.M(6,:)); title({sPlot.idStr '2 Sentences FFTs'});
figure; scatter3(mPhi(:,1), mPhi(:,2), mPhi(:,3), sPlot.pointSize, 1:size(sCsf.M,1), 'Fill'); colorbar;
sPlot.lineNums = cellfun(@num2str, num2cell(1:size(sCsf.M,1)), 'UniformOutput', false);
text(mPhi(:,1), mPhi(:,2), mPhi(:,3),sPlot.lineNums,'FontSize',sPlot.fontSize,'Color','red'); 
title({sPlot.idStr 'First 3 Eigen Vectors'} );
figure; scatter3(mPhi(:,2), mPhi(:,3), mPhi(:,4), sPlot.pointSize, 1:size(sCsf.M,1), 'Fill'); colorbar;
text(mPhi(:,2), mPhi(:,3), mPhi(:,4),sPlot.lineNums,'FontSize',sPlot.fontSize,'Color','red'); 
title({sPlot.idStr '2,3,4 Eigen Vectors'} );
% figure; plot(mPhi(:,1));
% figure; plot(mPhi(:,2));
% figure; plot(mPhi(:,3));

clear('tmp','i');