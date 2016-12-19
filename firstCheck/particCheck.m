
clearvars -except sRawXlDat ; close all;
addpath('../func');

%% sim parameters
extractData             = 1;
sXl.participantID       = '324.1';
featureName             = 'GazePointY' ;  % options are: 0:PupilRight, 1:PupilLeft, 2:GazePointX, 3:GazePointY .
sCsf.eigVecNum          = 4;
%sCsf.type               = 'pca';
sCsf.type               = 'diffusion maps';
sCsf.sDifMap.gausFactor = 8;
sCsf.minVecLength       = 50;

%% add participant data folder
sXl.dataFolderPath           = '../Data/Ordered/';
sXl.participantFolderPath    = strcat(sXl.dataFolderPath, sXl.participantID);
sXl.participantFilePathNoExt = strcat(sXl.participantFolderPath,'/',sXl.participantID);
sXl.xlFileName               = strcat(sXl.participantFilePathNoExt,'.xlsx');
addpath(sXl.participantFolderPath);
if extractData == 1
[sRawXlDat.num,sRawXlDat.text,sRawXlDat.raw]  = xlsread(sXl.xlFileName);
end;
%% relevent excel columns
sXl.sCol.PupilLeft  = 72; sXl.sCol.PupilRight = 73;
sXl.sCol.GazePointX = 57; sXl.sCol.GazePointY = 58; 
sXl.sCol.ValidLeft  = 74; sXl.sCol.ValidRight = 75;
sXl.sCol.EventType  = 37; sXl.sCol.EventName  = 38;
sXl.LineStart     = 'InstructionStart'; 
sXl.LineEnd       = 'InstructionEnd';
%% get idxs of lines - only valid data
sXlDat.mNum   = sRawXlDat.num;
sXlDat.mTxt   = sRawXlDat.text(2:end,:);
sXlDat.vValidEyesIdx = find((sXlDat.mNum(:,sXl.sCol.ValidRight) == 0 & sXlDat.mNum(:,sXl.sCol.ValidLeft) == 0) | ...
                        strcmp(sXlDat.mTxt(:,sXl.sCol.EventType),sXl.LineStart) | ...
                        strcmp(sXlDat.mTxt(:,sXl.sCol.EventType),sXl.LineEnd) );
sXlDat.mNumValid  = sXlDat.mNum(sXlDat.vValidEyesIdx,:);
sXlDat.mTxtValid  = sXlDat.mTxt(sXlDat.vValidEyesIdx,:);
sXlDat.vStartInd  = find(strcmp(sXlDat.mTxtValid(:,sXl.sCol.EventType),sXl.LineStart))+1;
sXlDat.vStartInd  = sXlDat.vStartInd(2:end);
sXlDat.vEndInd    = find(strcmp(sXlDat.mTxtValid(:,sXl.sCol.EventType),sXl.LineEnd))-1;
sXlDat.vEndInd    = sXlDat.vEndInd(2:end);
if length(sXlDat.vEndInd) ~= length(sXlDat.vStartInd)
    error('StartIdx and EndIdx must be the same size, check Xl file'); end
sXlDat.numOfSamp      = length(sXlDat.vEndInd);
%% get feature data
sXlDat.vFeatureIdx    = getfield(sXl.sCol, featureName);
sXlDat.vFeatureDat    = sXlDat.mNumValid(:,sXlDat.vFeatureIdx);
[sDat.M,sDat.minLength,sDat.maxLength] = split2vecs(sXlDat.vFeatureDat,sXlDat.vStartInd,sXlDat.vEndInd);
%% spectrum of data
if sDat.minLength < sCsf.minVecLength;
    sSpec.minLength = sCsf.minVecLength;
else
    sSpec.minLength = sDat.minLength;
end
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
sCsf.M            = sSpec.M;
if strcmp(sCsf.type,'diffusion maps')
    [mPhi, mLam] = gausDiffusinMaps(sCsf.M,'median',sCsf.sDifMap.gausFactor,sCsf.eigVecNum);
elseif strcmp(sCsf.type,'pca')
    [mPhi, mLam] = pcaWrap(sCsf.M,1,sCsf.eigVecNum);
end

%%
close all;
sPlot.idStr = [sXl.participantID '  ' featureName];
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