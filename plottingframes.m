%% Plotting for EMBC 
clc
clear all
TestFolders=[ "jan12"];

for iTest=1:length(TestFolders)
    TestFiles{iTest}=sprintf("%s_ana",TestFolders{iTest});
end

AnaStruct=sprintf("%s_ana",TestFolders);
S = load_test(TestFolders,TestFiles);
%% 
TestStruct=sprintf("%s_test",TestFolders);
AmpGain=S.(TestStruct).ExpPar.AmpGain;
Lbl='Occ';
ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
PlotTrial=[7 7 ];
PlotFrame=floor([ 489 491]);
ExpLabel=ExpTable{1};
DataLabels=S.(AnaStruct).AnaPar.DataLabels;
BlankLength=S.(AnaStruct).AnaPar.BlankLength;
fs=S.(TestStruct).ExpPar.fs;
FiltLabel="GS";
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
    ylim([1.2*min(BPFilt_EMG) 1.2*max(BPFilt_EMG)])
    Frames=PlotFrame(1):PlotFrame(2);

    for iFrame =1:length(Frames)
        x=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time')(BegofFrames(Frames(iFrame)));
        plot([ x x ]+1/fs,[-10 10],'--k','LineWidth',1);
        plot([ x x ]+BlankLength/fs,[-10 10],'--k','LineWidth',1);

    end
    
    grid on 
    subplot(2,1,2)
    MWave=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).mWaveswithDropped(:,PlotFrame(1):PlotFrame(2));
    vEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGwithDropped(:,PlotFrame(1):PlotFrame(2));
    plot(Time,[MWave(:,1); MWave(:,2); MWave(:,3)],'k','LineWidth',2)
    hold on
    plot(Time,[vEMG(:,1); vEMG(:,2); vEMG(:,3)],'r','LineWidth',2)
    grid on
    ttl=sprintf(" %s Estimated vEMG and MWaves",FiltLabel);
    title(ttl);
    xlabel("Time (s)")
    ylabel("EMG")
    ylim([1.2*min(BPFilt_EMG) 1.2*max(BPFilt_EMG)])
    lgd=legend({'MWaves','vEMG'});

end

%%
TestStruct=sprintf("%s_test",TestFolders);
AmpGain=S.(TestStruct).ExpPar.AmpGain;
Lbl='Occ';
ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
PlotTrial=[9 9 ];
PlotFrame=floor([ 487 489]);
ExpLabel=ExpTable{1};
DataLabels=S.(AnaStruct).AnaPar.DataLabels;
BlankLength=S.(AnaStruct).AnaPar.BlankLength;
fs=S.(TestStruct).ExpPar.fs;
FiltLabel="GS";
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

    figure(iTrial)
    subplot(2,1,1)
    plot(Time,EMG/AmpGain,'k','LineWidth',2)
    hold on
    plot(Time,BlankEMG,'r','LineWidth',2)
    lgd=legend({'Before Blanking and BP Filtering', 'After Blanking and BP Filtering'});
    lgd.AutoUpdate=false;
%     ttl=sprintf('BP Filtered %s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{1},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2));
%     title(ttl);
    title("An Example Dropped Frame")
    xlabel("Time (s)")
    ylabel("EMG")
    ylim([1.2*min(BPFilt_EMG) 1.2*max(BPFilt_EMG)])
    Frames=PlotFrame(1):PlotFrame(2);

    for iFrame =1:length(Frames)
        x=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time')(BegofFrames(Frames(iFrame)));
        plot([ x x ]+1/fs,[-10 10],'--k','LineWidth',1);
%         plot([ x x ]+BlankLength/fs,[-10 10],'--k','LineWidth',1);

    end

    grid on 
    text(13.65,4*10^-4,'Dropped Frame')

    S.(TestStruct).(ExpLabel).TrialsPW
    DataLabels=S.(AnaStruct).AnaPar.DataLabels;
    DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
    DroppedEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedEMG;
    clear lgd

    for iFrame=1:length(DroppedFrames)
        figure(iTrial)
        subplot(2,1,2)
        plot(DroppedEMG(:,iFrame),'LineWidth',2)
        hold on
        lgd(iFrame)=sprintf("Frame %d",DroppedFrames(iFrame));
        legend(lgd,'Location','NorthWest');
%         ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
%         title(ttl);
        title('EMG Signals of Dropped Frames')
        xlabel('Samples')
        ylabel('BP Filtered EMG')
    end
    grid on
end


%% Experiment figures 
close all

PlotTime=[0 10];
exp_lbl='RC';
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
stim_freq=S.(TestStruct).ExpPar.freq_list(1);
fs=S.(TestStruct).ExpPar.fs;
st=1/fs;
for iTest=1:1
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;

    PlotFrame=[PlotTime(1)*stim_freq PlotTime(2)*stim_freq];
    PlotInd=[PlotFrame(1):PlotFrame(2)];
    FiltLabel="Unfilt";
    clear lgd
    for iPW=5:5
        
        Ind_PW(:,iPW)=PWofTrials==PWPoints(iPW);
        IndTrials(:,iPW)=find(Ind_PW(:,iPW)==1);
        
        for iTrial=1:length(IndTrials(:,iPW))
            TrialLabel=sprintf("Trial_%d",IndTrials(iTrial,iPW));
        
            F=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Force');
            T=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time');
            F= F(T>PlotTime(1) & T<PlotTime(2));
            T=T(T>PlotTime(1) & T<PlotTime(2));
            figure(iPW)
            subplot(2,1,1)
            plot(T,F,'k','LineWidth',2)
            hold on
            lgd{iTrial}=sprintf("Trial %d",IndTrials(iTrial,iPW));
        end
        plot(PWProfile(1,2:5),PWProfile(2,2:5)*2,'r','LineWidth',2)

        grid on
        title('Stimulated Force during Recruitment Curve Trials')
        xlabel('Time(s)')
        ylabel('Force(N)')
        lgd{end+1}="PW (a.u.)";
        legend(lgd)
        xlim([0 10])
    end
end



PlotTime=[0 15];
vMVC = 1;
sMVC = 2;
exp_lbl='Occ';
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
VoliLevels=[10 20 30 40]; VoliLevel=VoliLevels(vMVC);
StimLevels=[ 0 10 20 30]; StimLevel=StimLevels(sMVC);

stim_freq=S.(TestStruct).ExpPar.freq_list(1);
fs=S.(TestStruct).ExpPar.fs;
st=1/fs;

for iTest=1:1
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
    PlotFrame=[PlotTime(1)*stim_freq PlotTime(2)*stim_freq];
    PlotInd=[PlotFrame(1):PlotFrame(2)];
    FiltLabel="Unfilt";
    IndTrials= find_trialnum(VoliLevel, StimLevel, RepTableMat);

    clear lgd

    for iTrial=1:length(IndTrials)
        TrialLabel=sprintf("Trial_%d",IndTrials(iTrial));

        Target=S.(TestStruct).(ExpLabel).TargetProfile2;
        F=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Force');
        PW=S.(TestStruct).(ExpLabel).(TrialLabel).data.('PW');
        T=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time');
        F= F(T>PlotTime(1) & T<PlotTime(2));
        PW=PW(T>PlotTime(1) & T<PlotTime(2));

        T=T(T>PlotTime(1) & T<PlotTime(2));
        figure(iPW)
        subplot(2,1,2)
        plot(T,F*0.95,'k','LineWidth',2)
        hold on
        lgd{iTrial}=sprintf("Trial %d",IndTrials(iTrial));
    end
    
    plot(Target(1,1:4),Target(2,1:4)*RepTableMat(IndTrials(iTrial),1),'r','LineWidth',2)
    lgd{end+1}="Target Line";
    legend(lgd)
    grid on
    title('Force Measurement Superimposed to Visual Target Line')
    xlabel('Time(s)')
    ylabel('Force(N)')
    
end
%% statistical tests results 





















