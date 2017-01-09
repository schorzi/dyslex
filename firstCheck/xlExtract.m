function [ sXl, sXlDat ] = xlExtract( datFoldPath, participantID )

mp3Fs = 44100;
sXl.participantIDnum = str2double(participantID); 

vParticipantsIDs = [252,293,295,298,323,324,326,327,328];
vSixLineStart    = [37337 22586 37000 0 0 0 0 42629 0 ];   % in micro seconds %252=19700 295=17221

%% add participant data folder
sXl.participantID            = participantID;
sXl.dataFolderPath           = datFoldPath;
sXl.participantFolderPath    = strcat(sXl.dataFolderPath, sXl.participantID,'/');
sXl.participantFilePathNoExt = strcat(sXl.participantFolderPath,sXl.participantID);
sXl.xlFileName               = strcat(sXl.participantFilePathNoExt,'.xlsx');
sXl.linesIdxFile             = strcat(sXl.participantFolderPath,'linesIdx.mat');

sXl.sixLineStart             = vSixLineStart(vParticipantsIDs==sXl.participantIDnum);
%% extract data from excel
[sXlDat.num,sXlDat.text,sXlDat.raw]  = xlsread(sXl.xlFileName);

load(sXl.linesIdxFile);  % variable linesIdx
sXlDat.mlinesTimes = linesIdx * (1000/mp3Fs); 
sXlDat.numOfLines = size(sXlDat.mlinesTimes,1);

%% relevent excel columns
sXl.sCol.('PupilLeft')  = 72; sXl.sCol.('PupilRight') = 73;
sXl.sCol.('GazePointX') = 50; sXl.sCol.('GazePointY') = 51; 
sXl.sCol.ValidLeft  = 74; sXl.sCol.ValidRight = 75;
sXl.sCol.EventType  = 37; sXl.sCol.EventName  = 38;
sXl.sCol.time       = 18; % in num mtx
sXl.LineStart       = 'InstructionStart'; 
sXl.LineEnd         = 'InstructionEnd';

%% get time
sXlDat.vTime            = sXlDat.num(:,sXl.sCol.time);    % time in milliseconds
sXlDat.vTime(isnan(sXlDat.vTime))=[];                      % remove empty
sXlDat.vTime          = (sXlDat.vTime - sXlDat.vTime(1));  % align start zero

sXl.xl2mp3delay       = sXl.sixLineStart - sXlDat.mlinesTimes(6,1);    % set mp3 to excel delay
sXlDat.mLinesTimeSync = sXlDat.mlinesTimes + sXl.xl2mp3delay;
% get idxs of closest values to time
sXlDat.edges          = [-Inf, mean([sXlDat.vTime(2:end)'; sXlDat.vTime(1:end-1)']), +Inf];
sXlDat.mLinesIdx      = discretize(sXlDat.mLinesTimeSync, sXlDat.edges);
sXlDat.mNum           = sXlDat.num;
sXlDat.mTxt           = sXlDat.text(2:end,:);

%% extract all features data : PupilRight, PupilLeft, GazePointX, GazePointY
sXlDat.('PupilRight') = sXlDat.mNum(:,sXl.sCol.PupilRight);
sXlDat.('PupilLeft')  = sXlDat.mNum(:,sXl.sCol.PupilLeft);
sXlDat.('GazePointX') = sXlDat.mNum(:,sXl.sCol.GazePointX);
sXlDat.('GazePointY') = sXlDat.mNum(:,sXl.sCol.GazePointY);
end

