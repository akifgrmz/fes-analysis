%% Plotting frames 
clc
clear all
TestFolders=[ "jan12"];

for iTest=1:length(TestFolders)
    TestFiles{iTest}=sprintf("%s_ana",TestFolders{iTest});
end

AnaStruct=sprintf("%s_ana",TestFolders);
S = load_test(TestFolders,TestFiles);
%%
TestStruct=sprintf("%s_test",TestFolders)
AmpGain=S.(TestStruct).ExpPar.AmpGain;
Lbl='Occ';
ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
PlotTrial=[7 7];
PlotFrame=floor([ 487 489]);
ExpLabel=ExpTable{1};
DataLabels=S.(AnaStruct).AnaPar.DataLabels;
BlankLength=S.(AnaStruct).AnaPar.BlankLength;
fs=S.(TestStruct).ExpPar.fs;
for iTrial=PlotTrial(1):PlotTrial(2)
    TrialLabel=sprintf('Trial_%d',iTrial);
    
    FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;
    BegofFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
    Ind= [BegofFrames(PlotFrame(1)) BegofFrames(PlotFrame(2)+1)-1];
    
    BPFilt_EMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BPFilt_EMG');
    BlankEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BlankedEMG');
    Trig=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Trigger');
    EMG=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('EMG');
    Time=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Time');

    figure
    subplot(2,1,1)
    plot(Time,Trig/2000,'b','LineWidth',2)
    hold on
    plot(Time,EMG/AmpGain,'k','LineWidth',2)
    plot(Time,BlankEMG,'r','LineWidth',2)
    lgd=legend({'Trigger (a.u.)','Before Blanking and BP Filtering', 'After Blanking and BP Filtering'});
    lgd.AutoUpdate=false;
    ttl=sprintf('BP Filtered %s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{1},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2));
    title(ttl);
    xlabel("Time (s)")
    ylabel("EMG")
    ylim([1.5*min(BPFilt_EMG) 1.5*max(BPFilt_EMG)])
    
    Frames=PlotFrame(1):PlotFrame(2);

    for iFrame =1:length(Frames)
        x=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time')(BegofFrames(Frames(iFrame)));
        plot([ x x ],[-10 10],'--k','LineWidth',1);
        plot([ x x ]+BlankLength/fs,[-10 10],'--k','LineWidth',1);

    end
    grid on 
    subplot(2,1,2)

end
