%% This script is for tidying the data after an experiment
%Run only one time
% This will reorganize the data for the following analyses 
clear all
%
% FolderNames='dec5'; 
% only one file at a time 
% FolderNames=["nov8", "nov28_2","dec5"];  % Foldername to be loaded
FolderNames=["dec5","nov28_2","nov27","nov8","nov7"]; 
iFolder=5;
FolderName=FolderNames(iFolder);
FileName='Expsave';     % File name
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

RCTurnOffTimes= table([ 10 ],[10],[ 10],[8],[10],...
    'VariableNames',["dec5","nov28_2","nov27","nov8","nov7"]);
S.(ExpStruct).RCCurveTrials.TurnOffTime=table2array(RCTurnOffTimes(:,FolderName));

mrg=1.8;
T=S.(ExpStruct).RCCurveTrials.Trial_1.data(1000:2000,5);
sample_t=abs((T(1)-T(end))/(2000-1000));
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


%% First plotting program
%File loading section

clear all

iTestPlot=1;  % pick a test to plot 
iExp=2;    % pick an exp to plot 
FolderNames={'nov8','dec5'};  %% Folders to be loaded 
FileNames={'nov8_test','dec5_test'};  %% Files to be loaded 

M = load_exp(FolderNames,FileNames);
Fields = fieldnames(M);
ExpStruct=Fields{iTestPlot};
ExpLabels=M.(ExpStruct).ExpLabels;
ExpLabel=ExpLabels{iExp};
DtInd= M.(ExpStruct).ExpPar.DataInd;
TimeRange=[4 22];  % in seconds
PlotRange=[ 10 10 ]; 

iForce=table2array(DtInd(:,"Force"));
iTrigger=table2array(DtInd(:,"Trigger"));
iEMG=table2array(DtInd(:,"EMG"));
iTime=table2array(DtInd(:,"Time"));
iPW=table2array(DtInd(:,"PW"));

%
for iTrial=PlotRange(1):PlotRange(2)
    TrialLabel=sprintf('Trial_%d',iTrial);
    
    Time=M.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iTime);
    TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
    
    EMG=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iEMG);
    Trigger=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iTrigger);
    Force=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iForce);
    PW=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iPW);
    figure(iTrial)
    subplot(2,1,1)
    plot(Time(TimeInd),Trigger/5,'r','LineWidth',1)
    hold
    plot(Time(TimeInd),EMG,'b','LineWidth',2)
    legend({'Trigger(a.u.)','Raw EMG (a.u.)'})
    title('Trigger and EMG Signal')
    xlabel('Time (s)')
    
    subplot(2,1,2)
    plot(Time(TimeInd),PW/20,'r','LineWidth',1)
    hold
    plot(Time(TimeInd),Force,'b','LineWidth',2)
    legend({'PW (a.u.)','Force(N)'})
    title('Pulse Width and Force Signal')
%     set(gca,'XTick',[5 8 11 12 15 18 21])
    xlabel('Time (s)')

%     legend({'Pre-BPFiltered', 'Post-BPFiltered'})
%     ttl=sprintf('%dth order BP filtered %s at %s, TrialNum: %d',BPOrder,DataLabels{iEMG},ExpLabel,iTrial);
%     title(ttl);
%     xlabel(DataLabels{iTime})
%     ylabel(DataLabels{iEMG})

end
M.(ExpStruct).OccTrials.RepTableMat

%% Checking the recruitment curves
% NOT READY
iTestPlot=2;  % pick a test to plot 
FolderNames={'nov28_2','dec5'};  %% Folders to be loaded 
iExp=2;
TimeRange=[4 22];  % in seconds
PlotRange=[ 10   15 ]; 
iForce=3;
iTrigger=2;
iEMG=1;
iTime=5;
iPW=4;

for iFile=1:length(FolderNames)
    ExpStruct=sprintf('%s_test',char(FolderNames{iFile}));
    str=sprintf('%s/%s',char(FolderNames{iFile}),ExpStruct);
    load (str);
    fields = fieldnames(S);
    M.(fields{1})=S.(fields{1});
end

% works only with dec5 for now 
ExpStruct=sprintf('%s_test',char(FolderNames{iTestPlot}));
ExpLabels=M.(ExpStruct).ExpLabels;
y=M.(ExpStruct).(ExpLabels{iTestPlot}).MeanForce;
RedoTrials=M.(ExpStruct).(ExpLabels{iTestPlot}).RedoTrials;
for iRedo=1:length(RedoTrials)
    RedoLabel=sprintf('RedoTrial_%d',iRedo);
    y(iRedo)=M.(ExpStruct).(ExpLabels{iTestPlot}).(RedoLabel).AvgForce;
end

%     gompertz



function Vars = gompertz(a, PW)

Vars = a(1)*exp(a(2)*exp(a(3)*PW));

end

