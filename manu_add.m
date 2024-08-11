%% manual additions 

clear all
%TestFolders=["feb28_24" ];
%TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25" "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=["feb28_24" "feb29_24" "mar18_24" "mar20_24"];
TestFolders=["jan7" "jan11" "jan12"   ];
% TestFolders=["jun20_24" ]

for iTest=1:length(TestFolders)
    TestFolder=TestFolders(iTest);
    S=tidy_data(TestFolder);
    TestFile=sprintf("%s_test",TestFolder);
    str=sprintf('%s/%s',TestFolder,TestFile);
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
end

%% 1- Jun20_24
% Load

FolderName=["jun20_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
% Additions 

S.(TestFile).ExpPar.DataInd=["EMG","Trigger","Force","PW","Time","Torque","UDP","Hand"];
S.(TestFile).ExpPar.ExpLabels=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials","RCRampTrials" ];
S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
S.(TestFile).ExpPar.EffortType="Hand";
S.(TestFile).OccTrials.StimRange=[5 11];
S.(TestFile).OccTrials.StimConstantRange=[6 11];
S.(TestFile).OccTrials.StimProfile=11; % Stim turn off time 
% Save

str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile)
str=sprintf("%s.mat is created",TestFile);
disp(str)

%% 2- jul9_24
% Load

FolderName=["jul9_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
% Additions 

S.(TestFile).ExpPar.DataInd=["EMG","Trigger","Force","PW","Time","Torque","UDP","Hand"];
S.(TestFile).ExpPar.ExpLabels=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials","RCRampTrials" ];
S.(TestFile).ExpPar.EffortType="Hand";
S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
S.(TestFile).OccTrials.StimRange=[5 11];
S.(TestFile).OccTrials.StimConstantRange=[6 11];
S.(TestFile).OccTrials.StimProfile=11; % Stim turn off time 
% Save

str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile) 
str=sprintf("%s.mat is created",TestFile);
disp(str)

%% 3- Jul21_24
% Load

FolderName=["jul21_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
% Additions 

S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
% Save

str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile)
str=sprintf("%s.mat is created",TestFile);
disp(str)

%% 4- Jul31_24
% Load
FolderName=["jul31_24" ];
TestFile=sprintf("%s_test",FolderName);
S=tidy_data(FolderName,false);
% Additions 

S.(TestFile).ExpPar.DataIndTable=array2table(1:length(S.(TestFile).ExpPar.DataInd),'VariableNames',S.(TestFile).ExpPar.DataInd);
S.(TestFile).OccTrials.StimProfile=14; % Stim turn off time 
S.(TestFile).OccTrials.StimRange=[3 14];
S.(TestFile).OccTrials.StimConstantRange=[7 14];
% Save

str=sprintf('%s/%s',FolderName,TestFile);
save(str,'-struct','S',TestFile)
str=sprintf("%s.mat is created",TestFile);
disp(str)