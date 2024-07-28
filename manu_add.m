%% manual additions 

%% 1- Jun20_24

FolderName=["jun20_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);

S.(TestFile).ExpPar.DataInd=["EMG","Trigger","Force","PW","Time","Torque","UDP","Hand"];
S.(TestFile).ExpPar.ExpLabels=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials","RCRampTrials" ];
S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
S.(TestFile).ExpPar.EffortType="Hand";
S.(TestFile).OccTrials.StimRange=[5 11];
S.(TestFile).OccTrials.StimProfile=11; % Stim turn off time 

str=sprintf('%s/%s',FolderName,TestFile)
save(str,'-struct','S',TestFile)
TestFile
str=sprintf("%s.mat is created",TestFile);
disp(str)

%% 2- jul9_24

FolderName=["jul9_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
S.(TestFile).ExpPar.DataInd=["EMG","Trigger","Force","PW","Time","Torque","UDP","Hand"];
S.(TestFile).ExpPar.ExpLabels=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials","RCRampTrials" ];
S.(TestFile).ExpPar.EffortType="Hand";
S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
% S.(TestFile).OccTrials.StimRange=[5 11];
% S.(TestFile).OccTrials.StimProfile=11; % Stim turn off time 

str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile) 
str=sprintf("%s.mat is created",TestFile);
disp(str)

%% 3- Jul21_24

FolderName=["jul21_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile)
str=sprintf("%s.mat is created",TestFile);
disp(str)

