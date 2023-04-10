%% First plotting program
%File loading section
clear all
close all
iTestPlot=1;  % pick a test to plot 
FolderNames={'jan11'};  %% Folders to be loaded 
FileNames={'jan11_test'};  %% Files to be loaded 


TimeRange=[4 22];  % in seconds
PlotRange=[ 3 3 ];  % Trial
iExp=5;    % pick an exp to plot 

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
    legend({'Trigger(a.u.)','Raw EMG (V)'})
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
