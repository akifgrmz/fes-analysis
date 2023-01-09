function tidy_data(FolderNames)

%% This script is for tidying the data after an experiment
% Enter a string of the folder name where raw data is located. raw file
% name must be "expsave.mat"
% This will reorganize the data for the following analyses 
% Only one file at a time 

FolderNames=["jan7"]; 
% FolderNames=["nov8", "nov28_2","dec5"];  % Foldername to be loaded
% FolderNames=["dec5","nov28_2","nov27","nov8","nov7"]; 
iFolder=1;
FolderName=FolderNames(iFolder);
FileName='expsave';     % File name
ExpStruct=sprintf('%s_test',FolderName);

str=sprintf('%s/%s.mat',FolderName,FileName);
temp=load (str);
S.(ExpStruct).ExpPar.Amp=temp.handles.Amp;
S.(ExpStruct).ExpPar.MaxPW=temp.handles.MaxPW;
S.(ExpStruct).ExpPar.PW=temp.handles.PW;
S.(ExpStruct).ExpPar.FreqList=temp.handles.freq_list;
S.(ExpStruct).ExpPar.CalibMatrix=temp.handles.CalibMatrix;
S.(ExpStruct).ExpPar.FolderName=temp.handles.FolderName;
S.(ExpStruct).ExpPar.sample_t=temp.handles.sample_t;
S.(ExpStruct).ExpPar.fs=1/temp.handles.sample_t;
S.(ExpStruct).ExpPar.FolderName=temp.handles.FolderName;
S.(ExpStruct).ExpPar.freq_list=temp.handles.freq_list;
S.(ExpStruct).ExpPar.CalibMatrix=temp.handles.CalibMatrix;
S.(ExpStruct).ExpPar.DataIndTable=temp.handles.DataIndTable;


S.(ExpStruct).MVCTrials=temp.handles.MVCTrials;
S.(ExpStruct).FatigueTrials=temp.handles.FatigueTrials;
S.(ExpStruct).CustomTrials=temp.handles.CustomTrials;
S.(ExpStruct).OccTrials=temp.handles.ExpTrials;
S.(ExpStruct).RCCurveTrials=temp.handles.RCCurveTrials;


%
% Labeling
%

DataInd={'EMG','Trigger','Force','PW','Time'};
ExpLabels={'MVCTrials','RCCurveTrials','CustomTrials','OccTrials','FatigueTrials' };
S.(ExpStruct).ExpPar.ExpLabels=ExpLabels;
S.(ExpStruct).ExpLabels=ExpLabels;
S.(ExpStruct).ExpPar.DataInd=table(1,2,3,4,5,'VariableNames',DataInd);


% 
% Num of Trials might not be needed anymore
%
for iExp=1:length(ExpLabels)
    ExpLabels=S.(ExpStruct).ExpLabels;
    ExpLabel=ExpLabels{iExp};
    
    NumofTrials=S.(ExpStruct).(ExpLabel).iTrial-1;
    S.(ExpStruct).(ExpLabel).NumofTrials=NumofTrials;
    c=1;
    TrialNums=[];
    S.(ExpStruct).(ExpLabel).RedoTrials=TrialNums;
    for iTrial=1:NumofTrials
        RedoLabel=sprintf('RedoTrial_%d',iTrial);
        
        if isfield(S.(ExpStruct).(ExpLabel), RedoLabel)
            TrialNums(c)=iTrial;
            c=c+1;
        end
        S.(ExpStruct).(ExpLabel).RedoTrials=TrialNums;
        S.(ExpStruct).(ExpLabel).NumofRedos=length(TrialNums);
    end
    S.(ExpStruct).ExpPar.NumofTrialsTable{iExp,1}=ExpLabel;
    S.(ExpStruct).ExpPar.NumofTrialsTable{iExp,2}=S.(ExpStruct).(ExpLabel).NumofTrials;
    
end



%
% RC turn off times and data length
%

% RCTurnOffTimes= table([ 10 ],[10],[ 10],[8],[10],...
%     'VariableNames',["dec5","nov28_2","nov27","nov8","nov7"]);

TurnOffInd=4;
S.(ExpStruct).RCCurveTrials.TurnOffInd=TurnOffInd;
S.(ExpStruct).RCCurveTrials.TurnOffTime=S.(ExpStruct).RCCurveTrials.PWProfile(1,TurnOffInd);

mrg=1.8;
% T=S.(ExpStruct).RCCurveTrials.Trial_1.data(1000:2000,5);
% sample_t=abs((T(1)-T(end))/(2000-1000));
sample_t=S.(ExpStruct).ExpPar.sample_t;
fs=1/sample_t;

T_turnoff=S.(ExpStruct).RCCurveTrials.TurnOffTime;
lgt=round(T_turnoff*fs+fs*mrg);


for iExp=2:2
    ExpLabels=S.(ExpStruct).ExpLabels;
    ExpLabel=ExpLabels{iExp};
    
    NumofTrials=S.(ExpStruct).(ExpLabel).NumofTrials;
    RedoTrials=S.(ExpStruct).(ExpLabel).RedoTrials;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        [lg_dt,~]=size(S.(ExpStruct).(ExpLabel).(TrialLabel).data);
        if  lg_dt>=lgt
            S.(ExpStruct).(ExpLabel).(TrialLabel).data(lgt:end,:)=[];
        else
            S.(ExpStruct).(ExpLabel).(TrialLabel).data(end+1:lgt,:)=0;
        end

    end 
    
    for iTrial=1:length(RedoTrials)
        RedoLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));
        
        [lg_dt,~]=size(S.(ExpStruct).(ExpLabel).(RedoLabel).data);
        if  lg_dt>=lgt
            S.(ExpStruct).(ExpLabel).(RedoLabel).data(lgt:end,:)=[];
        else
            S.(ExpStruct).(ExpLabel).(RedoLabel).data(end+1:lgt,:)=0;
        end    
    end
    
    S.(ExpStruct).ExpPar.lgt(iExp)=lgt;

end

S.(ExpStruct).ExpPar.fs_calc=fs;
S.(ExpStruct).ExpPar.sample_t_calc=sample_t;



%
%Fixing PW invalid values issue (-inf)
%

DtInd=S.(ExpStruct).ExpPar.DataInd;
iPW=table2array(DtInd(:,"PW"));
for iExp=1:length(ExpLabels)
    ExpLabels=S.(ExpStruct).ExpLabels;
    ExpLabel=ExpLabels{iExp};
    
     NumofTrials=S.(ExpStruct).(ExpLabel).NumofTrials;
    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        x=S.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iPW);
        x(x<0)=0;
        S.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iPW)=x;
    end
    
    RedoTrials=S.(ExpStruct).(ExpLabel).RedoTrials;

    for iTrial=1:length(RedoTrials)
        TrialLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));
        
        x=S.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iPW);
        x(x<0)=0;
        S.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iPW)=x;
    end
   
end


%
% Saving as struc
%

str=sprintf('%s/%s',FolderName,ExpStruct);
save(str,'-struct','S',ExpStruct)

str=sprintf("%s.mat is created",ExpStruct);
disp(str)
end
