
clearvars -except sRawXlDat ; close all;
addpath('../func');
addpath('../func/MinMaxFilter');

%% sim parameters
participantID       = '327';
%% add participant data folder
dataFolderPath           = '../../Data/Ordered/';
participantFolderPath    = strcat(dataFolderPath, participantID);
participantFilePathNoExt = strcat(participantFolderPath,'/',participantID);
mp3FileName              = strcat(participantFilePathNoExt,'.mp3');
addpath(participantFolderPath);
[stereo,Fs] = audioread(mp3FileName);
mp3         = stereo(:,1);
player      = audioplayer(mp3,Fs);
%% split lines
threshold            = 0.02;
mp3Filt    = minmaxfilt1(abs(mp3), 10000,'max','same');  
% mp3Filt2   = minmaxfilt1(mp3Filt, 2500,'min','same');
% mp3Filt3   = minmaxfilt1(mp3Filt2, 25000,'max','same');
mp3bin     = single(mp3Filt > threshold);
mp3binfilt = medfilt1(mp3bin,2000);
mp3Diff    = filter([1 -1],1,mp3binfilt);
indStartTmp  = find(mp3Diff > 0);
indEndTmp    = find(mp3Diff < 0 );
%%

close all;
figure; plot(abs(mp3));
figure; plot(mp3Filt);
figure; plot(mp3bin);
figure; plot(mp3(1:300000));
figure; plot(mp3Filt(1:300000));
figure; plot(mp3Filt2(1:300000));
figure; plot(mp3Filt3(1:300000));
figure; plot(mp3bin); ylim([-0.5 1.5]);
figure; plot(mp3binfilt(1:200000)); ylim([-0.5 1.5]);
figure; plot(mp3Diff(1:200000));
%%
indStartTmp2 = indStartTmp(1:2:end);
indEndTmp2 = indEndTmp(1:2:end);
indStart = indStartTmp2;
indStart = indStartTmp;
indEnd = indEndTmp2;
indStartBack = round(indStart - (0.2*Fs));
%% manual editing - 252
%{
addEndIdx = [63 82 85];
indEnd(addEndIdx) = indEnd(addEndIdx) + 0.2*Fs; % add 0.2 s to end 63, 82, 85.
indEnd(83) = indEnd(83) + 0.4*Fs;  % add 0.4 s to 83
indEnd(67) = indEnd(68);% combine 67-68
indEnd(73) = indEnd(74);% combine 73-74 , 79-80
indEnd(79) = indEnd(80);
delInd = [6,10,68,80,74];
indEnd(delInd) = []; indStartBack(delInd) = [];
indStartBack(61) = indStartBack(62)-200000;  % add idx
indEnd(61) = indEnd(62)-180000;

%% manual editing - 295

indEnd(16) = indEnd(17); % delete : 1,3,20
indEnd(28) = indEnd(29); % combine 16-17, 28-29, 40-41, 65-66, 68-69
indEnd(40) = indEnd(41);
indEnd(65) = indEnd(66);
indEnd(68) = indEnd(69);
delInd = [1 3 20 17 29 41 66 69];
indEnd(delInd) = []; indStartBack(delInd) = [];
indEnd = indEnd-10000;
%}
%% manual editing - 295
indStartBack(1)   = 1;
indStartBack(end) = [];
indStartBack(end) = []; indEnd(end) = [];

indEnd(37)        = indEnd(37)-50000; % combine 38 with end of 37
indStartBack(38)  = indEnd(37);

indEnd(51) = indEnd(52);% combine 51-52 
                       
indEnd(43:end+1)       = indEnd(42:end);             % split 41 to 2 sentences
indStartBack(43:end+1) = indStartBack(42:end);  
indEnd(42)             = indEnd(41);
indEnd(41)             = (indEnd(41)+indStartBack(41))/2;
indStartBack(42)       = indEnd(41);

indEnd(45:end+1)       = indEnd(44:end);             % split 42 to 2 sentences (now 43)
indStartBack(45:end+1) = indStartBack(44:end);  
indEnd(44)             = indEnd(43);
indEnd(43)             = (indEnd(43)+indStartBack(43))/2;
indStartBack(44)       = indEnd(43);

indEnd(62:end+1)       = indEnd(61:end);           % split 58 add some time after (now 60)
indStartBack(62:end+1) = indStartBack(61:end);
indEnd(61)             = indEnd(60) + 10000;
indEnd(60)             = (indEnd(60)+indStartBack(60))/2;
indStartBack(61)       = indEnd(60);

indEnd(70:end+1)       = indEnd(69:end);            % split 65 (now 68)
indStartBack(70:end+1) = indStartBack(69:end);  
indEnd(69)             = indEnd(68);
indEnd(68)             = ((indEnd(68)+indStartBack(68))/2) - 10000;
indStartBack(69)       = indEnd(68);                                  

indEnd(54)=[]; indStartBack(54)=[];           % remove 52 (now 54

indEnd = indEnd - 10000;
indEnd(41) = indEnd(41) + 8000;
indEnd(43) = indEnd(43) + 8000;
indEnd(44) = indEnd(44) + 8000;
indEnd(47) = indEnd(47) + 8000;
indEnd(49) = indEnd(49) + 8000; indEnd(60) = indEnd(60) + 15000;
indEnd(59) = indEnd(59) + 8000; indEnd(71) = indEnd(71) + 8000;
indEnd(61) = indEnd(61) + 8000; indEnd(67) = indEnd(67) + 8000;

%%  327 lines split
participantID            = '327';
dataFolderPath           = '../../Data/Ordered/';
filePath    = strcat(dataFolderPath,participantID,'/lines/lines.xlsx');
linesIdx  = xlsread(filePath,'A1:B80');
indStartBack = linesIdx(:,1);
indEnd =  linesIdx(:,2);
%%
path = strcat('../../Data/Ordered/',participantID);
for i=1:length(indStartBack)
    audiowrite(strcat(path,'/lines/',num2str(i),'.ogg'),...
                mp3(indStartBack(i):indEnd(i)),Fs); i
end
%%
linesIdx = [indStartBack indEnd];
save(strcat(path,'/linesIdx.mat'),'linesIdx');
