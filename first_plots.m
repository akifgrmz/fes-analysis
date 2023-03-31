%% First plotting program
%File loading section
clear all
close all
iTestPlot=1;  % pick a test to plot 
FolderNames={'mar16'};  %% Folders to be loaded 
FileNames={'mar16_test'};  %% Files to be loaded 


TimeRange=[4 22];  % in seconds
PlotRange=[ 3 3 ];  % Trial
iExp=4;    % pick an exp to plot 

M = load_test(FolderNames,FileNames);
Fields = fieldnames(M);
ExpStruct=Fields{iTestPlot};
ExpLabels=M.(ExpStruct).ExpPar.ExpLabels;
ExpLabel=ExpLabels{iExp};
DataInd= M.(ExpStruct).ExpPar.DataInd;

iForce=table2array(DataInd(:,"Force"));
iTrigger=table2array(DataInd(:,"Trigger"));
iEMG=table2array(DataInd(:,"EMG"));
iTime=table2array(DataInd(:,"Time"));
iPW=table2array(DataInd(:,"PW"));
TableInd=DataInd.Properties.VariableNames;

for iTrial=PlotRange(1):PlotRange(2)
    TrialLabel=sprintf('Trial_%d',iTrial);
    
    Time=M.(ExpStruct).(ExpLabel).(TrialLabel).data.(TableInd{iTime});
    TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
    
    EMG=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).(TableInd{iEMG});
    Trigger=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).(TableInd{iTrigger});
    Force=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).(TableInd{iForce});
    PW=M.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).(TableInd{iPW});
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
OccTable=M.(ExpStruct).OccTrials.RepTableMat;
OccTable=[[1:length(OccTable)]' OccTable]
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

%  gompertz



function Vars = gompertz(a, PW)

Vars = a(1)*exp(a(2)*exp(a(3)*PW));

end