%% Here Tidy up for Rstudio
% Are redos included here?
clear all
FolderNames={'dec5'};

for iFile=1:length(FolderNames)
    FileName=sprintf("%s_test",FolderNames{iFile}); 
    str=sprintf('%s/%s.mat',FolderNames{iFile},FileName);
    load (str);
    ExpLabels=S.(FileName{1}).ExpLabels;
    
    VarNames={'EMG','Trigger','Force','PW','Time','ExpNum','TrialNum'};
    sz = [1 length(VarNames)];
    VarTypes = {'double','double','double','double','double','string','string'};
    T = table('Size',sz,'VariableTypes',VarTypes,'VariableNames',VarNames);

    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels{iExp};
        NumofTrials=S.(FileName{1}).(ExpLabel).NumofTrials;
        %Here Fix the inconsistent data length
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            TrialData= S.(FileName{1}).(ExpLabel).(TrialLabel).data(:,1:5);

            Temp=array2table(TrialData,'VariableNames',VarNames(1:5));

            [row, ~] = size(TrialData);
            Exp_T=array2table( ExpLabel + strings(row,1),'VariableNames',VarNames(6));
            Trial_T=array2table( TrialLabel + strings(row,1),'VariableNames',VarNames(7));
            Temp(:,6)=Exp_T;
            Temp(:,7)=Trial_T;
            Temp = renamevars(Temp,["Var6","Var7"],[VarNames(6),VarNames(7)]);
            T=vertcat(T,Temp);
            clear Temp
        end
        
        RedoTrials=S.(FileName{1}).(ExpLabel).RedoTrials;
        if ~isempty(RedoTrials)
            
            for iRedo=1:length(RedoTrials)
                TrialLabel=sprintf('RedoTrial_%d',RedoTrials(iRedo));
                TrialData= S.(FileName{1}).(ExpLabel).(TrialLabel).data(:,1:5);

                Temp=array2table(TrialData,'VariableNames',VarNames(1:5));

                [row, ~] = size(TrialData);
                Exp_T=array2table( ExpLabel + strings(row,1),'VariableNames',VarNames(6));
                Trial_T=array2table( TrialLabel + strings(row,1),'VariableNames',VarNames(7));
                Temp(:,6)=Exp_T;
                Temp(:,7)=Trial_T;
                Temp = renamevars(Temp,["Var6","Var7"],[VarNames(6),VarNames(7)]);
                T=vertcat(T,Temp);
                clear Temp
            end
        end

    end
    T(1,:)=[];
    str=sprintf('%s/%sRtest.csv',FolderNames{iFile},FolderNames{iFile});
    writetable(T, str)
    height(T)
end
