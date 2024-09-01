function manu_add()
%% manual additions 

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
    
    %% 5- May31_24
    % Load
    
    FolderName=["may31_24" ];
    
    FileName='expsave';    
    ExpStruct=sprintf("%s_test",FolderName);
    AnaStruct=sprintf('%s_ana',FolderName);
    TestFile=sprintf("%s_test",FolderName);
    
    str=sprintf('%s/%s.mat',FolderName,FileName);
    temp=load (str);
    
    % Additions 
    
    temp.handles.ExpTrials.iTrial=16;
    str=sprintf('%s/%s',FolderName,'expsave');
    save(str,'-struct','temp','handles')
    
    S=tidy_data(FolderName,false);
    S.(TestFile).OccTrials.ListedNumofTrials=[16];
    
    S.(TestFile).OccTrials.StimProfile=11; % Stim turn off time 
    S.(TestFile).OccTrials.StimRange=[ 5 11];
    S.(TestFile).OccTrials.StimConstantRange=[6 11];
    S.(TestFile).OccTrials.dropped=true;
    
    % S.(TestFile).OccTrials.RepTableMat=S.(TestFile).OccTrials.RepTableMat(1:16,:);
    
    %Replace some trials
    for iTrial=1:3
        TrialLabel1=sprintf("Trial_%d",iTrial);
        TrialLabel2=sprintf("Trial_%d",iTrial+16);
    
        S.(TestFile).OccTrials.(TrialLabel1)=S.(TestFile).OccTrials.(TrialLabel2);
        S.(TestFile).OccTrials.(TrialLabel1).iTrial=iTrial+1;
    end
    
    % Save
    
    str=sprintf('%s/%s',FolderName,TestFile);
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
    
    %% 6- feb29_24
    
    FolderName=["feb29_24" ];
    TestFile=sprintf("%s_test",FolderName);
    S=tidy_data(FolderName,false);
    
    iTrial=7;
    TrialLabel1=sprintf("Trial_%d",iTrial);
    iEMG=S.(TestFile).ExpPar.DataIndTable.('EMG');
    Division=10;
    S.(TestFile).OccTrials.(TrialLabel1).data(:,iEMG)=S.(TestFile).OccTrials.(TrialLabel1).data(:,iEMG)/Division;
    
    str=sprintf('%s/%s',FolderName,TestFile);
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
    
    %% 7- aug19_24
    
    FolderName=["aug19_24"];
    TestFile=sprintf("%s_test",FolderName);
    S=tidy_data(FolderName,false);
    
    % Additions 
    S.(TestFile).ExpRuns(1)=true;
    S.(TestFile).MVCTrials.ExpRun=true;
    S.(TestFile).ExpRuns(6)=true;
    S.(TestFile).RCRampTrials.ExpRun=true;
    S.(TestFile).RCRampTrials.NumofTrials=1;
    S.(TestFile).RCRampTrials.iTrial=2;
    S.(TestFile).MVCTrials.NumofTrials=1;
    S.(TestFile).MVCTrials.iTrial=2;
    
    str=sprintf('%s/%s',FolderName,TestFile);
    
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
    
    %% 8- aug22_24
    
    FolderName=["aug22_24"];
    TestFile=sprintf("%s_test",FolderName);
    S=tidy_data(FolderName,false);
    
    % Additions 
    str=sprintf('%s/%s',FolderName,TestFile);
    
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
    %% 9- aug26_24
    
    FolderName=["aug26_24"];
    TestFile=sprintf("%s_test",FolderName);
    S=tidy_data(FolderName,false);
    
    % Additions 
    str=sprintf('%s/%s',FolderName,TestFile);
    
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
    
    %% 10- aug29_24
    
    FolderName=["aug29_24"];
    TestFile=sprintf("%s_test",FolderName);
    S=tidy_data(FolderName,false);
    
    % Additions 
    str=sprintf('%s/%s',FolderName,TestFile);
    
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)

end