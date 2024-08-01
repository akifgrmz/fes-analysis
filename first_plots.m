%% First plotting program
%File loading section
clear all
close all
iTestPlot=1;  % pick a test to plot 
FolderNames={'jan7'};  %% Folders to be loaded 
FileNames={'jan7_test'};  %% Files to be loaded 

% TestFiles(iTest)=sprintf("%s_test",TestFolders{iTest});
iExp=4;    % pick an exp to plot 

M = load_test(FolderNames,FileNames);
Fields = fieldnames(M);
ExpStruct=Fields{iTestPlot};
ExpLabels=M.(ExpStruct).ExpPar.ExpLabels;
%%
FolderNames=['mar20_24'];  %% Folders to be loaded 

PlotRange=[ 13 14 ];  % Trial
TimeRange=[4 22];  % in seconds
ExpStruct=sprintf('%s_test',FolderNames);
DataInd= S.(ExpStruct).ExpPar.DataIndTable;
iForce=DataInd.("Force");
iTrigger=(DataInd.("Trigger"));
iEMG=(DataInd.("EMG"));
iTime=(DataInd.("Time"));
iPW=(DataInd.("PW"));
TableInd=DataInd.Properties.VariableNames;

close all
for iTrial=PlotRange(1):PlotRange(2)
    TrialLabel=sprintf('Trial_%d',iTrial);
    ExpLabel=S.(ExpStruct).ExpPar.ExpLabels(4);

    Time=S.(ExpStruct).(ExpLabel).(TrialLabel).data(:,iTime);
    TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
    
    EMG=S.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iEMG);
    Trigger=S.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iTrigger);
    Force=S.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iForce);
    PW=S.(ExpStruct).(ExpLabel).(TrialLabel).data(TimeInd,iPW);
    figure(iTrial)
    subplot(2,1,1)
    plot(Time(TimeInd),Trigger/2,'r','LineWidth',1)
    hold
    plot(Time(TimeInd),EMG,'b','LineWidth',2)
    legend({'Trigger(a.u.)','Raw EMG (mV)'})
    title('Trigger and EMG Signal')
    xlabel('Time (s)')
    
    subplot(2,1,2)
    plot(Time(TimeInd),PW,'r','LineWidth',1)
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
OccTable=S.(ExpStruct).OccTrials.RepTableMat;
OccTable=[[1:length(OccTable)]' OccTable]
