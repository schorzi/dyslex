
clearvars -except sRawXlDat ; close all;
addpath('../func');

%% sim parameters

sCsf.featureName        = 'GazePointX' ;  % options are: 0:PupilRight, 1:PupilLeft, 2:GazePointX, 3:GazePointY .
sCsf.eigVecNum          = 4;
%Csf.type               = 'pca';
sCsf.type               = 'diffusion maps';
sCsf.sDifMap.gausFactor = 8;
sCsf.minVecLength       = 50;
%% relevent excel columns
sXl.sCol.PupilLeft  = 72; sXl.sCol.PupilRight = 73;
sXl.sCol.GazePointX = 57; sXl.sCol.GazePointY = 58; 
sXl.sCol.ValidLeft  = 74; sXl.sCol.ValidRight = 75;
sXl.sCol.EventType  = 37; sXl.sCol.EventName  = 38;
sXl.LineStart     = 'InstructionStart'; 
sXl.LineEnd       = 'InstructionEnd';
%%
vIds = dir('../Data/Ordered/');
sC.participantsID = {vIds(3:end).name};
specAll = [];
%%
for parNum=1:length(sC.participantsID)
    clearvars -except specAll sC sCsf parNum sXl
    %% add participant data folder
    parId = char(sC.participantsID(parNum));
    sC.dataFolderPath           = '../../Data/Ordered/';
    sC.participantFolderPath    = strcat(sC.dataFolderPath, parId);
    sC.participantFilePathNoExt = strcat(sC.participantFolderPath,'/',parId);
    sC.xlFileName               = strcat(sC.participantFilePathNoExt,'.xlsx');
    addpath(sC.participantFolderPath);
    [sRawXlDat.num,sRawXlDat.text,sRawXlDat.raw]  = xlsread(sC.xlFileName);
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
    sXlDat.vFeatureIdx    = getfield(sXl.sCol, sCsf.featureName);
    sXlDat.vFeatureDat    = sXlDat.mNumValid(:,sXlDat.vFeatureIdx);
    [sDat.M,sDat.minLength,sDat.maxLength] = split2vecs(sXlDat.vFeatureDat,sXlDat.vStartInd,sXlDat.vEndInd);
    %% spectrum of data
    sSpec.M = zeros(sXlDat.numOfSamp*sCsf.minVecLength,1);
    for i=1:sXlDat.numOfSamp
        tmp  = sDat.M(~isnan(sDat.M(:,i)),i);
        if length(tmp) < sCsf.minVecLength
            sSpec.M((i*sCsf.minVecLength)+1:(i+1)*sCsf.minVecLength)  = zeros(1,sCsf.minVecLength);
        else
            sSpec.M((i*sCsf.minVecLength)+1:(i+1)*sCsf.minVecLength)  = abs(fft(tmp,sCsf.minVecLength));
        end
    end
    
    specAll(parNum,:) = sSpec.M;
end
 %%
 validIdx = specAll(1,:)>0;
 for k=[2 3 4 8 9]
    validIdx = validIdx & specAll(k,:)>0;
 end
Idx = find(validIdx);
%%
 sCsf.M            = specAll([1 2 3 4 8 9],Idx);
if strcmp(sCsf.type,'diffusion maps')
    [mPhi, mLam] = gausDiffusinMaps(sCsf.M,'median',sCsf.sDifMap.gausFactor,sCsf.eigVecNum);
elseif strcmp(sCsf.type,'pca')
    [mPhi, mLam] = pcaWrap(sCsf.M,sCsf.eigVecNum,1);
end

%%
close all;
sPlot.idStr = [sCsf.featureName];
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