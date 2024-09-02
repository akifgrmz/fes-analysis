%% Manual fixes and additions

% Run once for each experiment
manu_add;


%% Tests and Parameters 
clc
clear all
% TestFolders=["mar20_24"];
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25" "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=[ "apr20" "oct11" "oct18"];
% TestFolders=["oct25"];
% TestFolders=[ "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=["jan7" "jan11" "jan12"  ];
% TestFolders=["jun20_24" "jul9_24" "jul21_24"  ]
% TestFolders=[ "jan7" "jan11" "jan12" "jun20_24" "jul9_24" "jul21_24" "jul31_24" ];
TestFolders=["feb28_24" "feb29_24" "mar18_24"  "mar20_24" ];
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24"];
TestFolders=[ "jan7" "jan11" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" "jan12" "aug22_24" "aug26_24"];
TestFolders=["aug29_24"];
  
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24" "aug29_24"];

NumofTests=length(TestFolders);
DroppedFrameFilts=strings(1,NumofTests)+"GS";
NoStimFilts=strings(1,NumofTests)+"Unfilt";
MAV_MAXMethods=strings(1,NumofTests)+"Real";

TauTestsForce=["jan7" "jan11" "jan12"]; % Average of 
TauTestsHand=["jan7" "jan11" "jan12"]; % Average of 

AvgOcclusionTests=["feb28_24"]; % Generalized result will not match the 
% indiv result for feb28_24 for this method because the slopes and averages will be slitly different
AvgOcclusionTests=["jan7" "jan11" "jan12";
                    "jun20_24" "jul9_24"	"jul21_24"]; 

% AvgOcclusionTests=["aug22_24" "aug26_24" "aug29_24" ];
AvgOcclusionTests=["jan7" "jan11" "jan12";
                     "aug22_24" "aug26_24" "aug29_24"]; 
ModelLog=false;
BlankTime=ones(1,NumofTests)*0.0035;
gs_orders=ones(1,NumofTests)*6;

S=ana_func(TestFolders,DroppedFrameFilts,NoStimFilts,MAV_MAXMethods,TauTestsForce,TauTestsHand,AvgOcclusionTests,ModelLog,BlankTime,gs_orders);
%% Save

%% Plots before analysis
clc
close all
Trials=[ 6 7 8]; 
TimeRange=[1 16];  % in seconds
% TestLabel=sprintf('%s_test',FolderNames);
ylims=5;
ExpNum=4;
for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    TestName=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpPar.ExpLabels(ExpNum);

    for iTrial=1:length(Trials)
        TrialLabel=sprintf('Trial_%d',Trials(iTrial));
        DataIndTable= S.(TestLabel).ExpPar.DataIndTable;

        Time=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,DataIndTable.("Time"));
        TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
        EffortType=S.(TestLabel).ExpPar.EffortType;
        EMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("EMG"));
        Trigger=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("Trigger"));
        Force=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.(EffortType));
        PW=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("PW"));

        figure(iTest)
        subplot(2,length(Trials),iTrial)
        plot(Time(TimeInd),Trigger/20,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),EMG,'b','LineWidth',2)
        legend({'Trigger(a.u.)','Raw EMG (mV)'})
        title(sprintf('Test:%s, Trial: %d', TestName, Trials(iTrial)))
        xlabel('Time (s)')
%         ylim([-ylims ylims])

        subplot(2,length(Trials),length(Trials)+iTrial)
        plot(Time(TimeInd),PW/200,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),Force,'b','LineWidth',2)
        legend({'PW (a.u.)','Meas. Effort'})
        title('Pulse Width and Force Signal')
        xlabel('Time (s)')
    end
OccTable=S.(TestLabel).OccTrials.RepTableMat;
OccTable=[[1:length(OccTable)]' OccTable]

end
%% Plots after filtering and preproccess 
clc
close all
Trials=[ 11 12 13];  % Trial
TimeRange=[1 12];  % in seconds
ylims=5;
ExpNum=4;
for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    TestName=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpPar.ExpLabels(ExpNum);
    EffortType=S.(TestLabel).ExpPar.EffortType;

    for iTrial=1:length(Trials)
        TrialLabel=sprintf('Trial_%d',Trials(iTrial));
        DataIndTable= S.(TestLabel).ExpPar.DataIndTable;

        TimeRange=S.(TestLabel).(ExpLabel).StimRange;
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
        TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);

        % RawEMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('EMG');
        BlankedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('BlankedEMG');

        Trigger=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('Trigger');
        EffortMeas=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).(EffortType);
        PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('PW');

        figure(iTest)
        subplot(length(Trials),1,iTrial)
        plot(Time(TimeInd),Trigger/2000,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),BlankedEMG,'b','LineWidth',2)
        legend({'Trigger(a.u.)','Raw EMG (mV)'})
        title(sprintf('Test:%s, Trial: %d', TestName, Trials(iTrial)))
        xlabel('Time (s)')
        % ylim([-1 1]*10^-3)
    end
OccTable=S.(TestLabel).OccTrials.RepTableMat;
OccTable=[[1:length(OccTable)]' OccTable]

end


%% Plotting, Debugging 
close all
Lbl='Occ';
PlotTrial=[ 10 10];
PlotTime=[ 12 13];
% PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);
PlotFrame=floor([ 478 482]);

TeststoPlot=TestFolders();

for iTest=1:length(TeststoPlot)
    TestLabel=sprintf("%s_test",TeststoPlot(iTest));
    AnaLabel=sprintf("%s_ana",TeststoPlot(iTest));
    DataLabels=S.(AnaLabel).AnaPar.DataLabels;
    TestName=S.(AnaLabel).AnaPar.TestName;
    ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).(Lbl));

    % PlotFrame=floor([ 300 305]); % first frame is skipped
    % ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    % ExpstoAna=ExpLabels([2,3,4,5]);
    FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;

    clear EMGFrames MWaveFrames vEMG MWave MWaveFrames_filtdropped vEMGFrames_filtdropped vEMG_filtdropped MWave_filtdropped
    % vEMGFrames_filtdropped=zeros(FrameLength,PlotFrame(2)-PlotFrame(1)+1,length(FiltLabels));
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameLength=S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength;

        for iFilt=1:length(FiltLabels)
            FiltLabel=FiltLabels{iFilt};

            EMGFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));

            vEMG(:,iFilt)=reshape(EMGFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave(:,iFilt)=reshape(MWaveFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

            vEMGFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));

            vEMG_filtdropped(:,iFilt)=reshape(vEMGFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave_filtdropped(:,iFilt)=reshape(MWaveFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

            % vEMGFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            % MWaveFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            % 
            % vEMG_withdropped(:,iFilt)=reshape(vEMGFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            % MWave_withdropped(:,iFilt)=reshape(MWaveFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        end

        TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
        TriggerFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TriggerFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));

        Time=reshape(TimeFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        Trig=reshape(TriggerFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{2},TeststoPlot(iTest)),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,Trig/10000,'b','LineWidth',2)
        hold on
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        
        plot(Time,vEMG(:,2),'r','LineWidth',2)
        legend({'Trigger(a.u)', 'Unfiltered', 'Comb vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d'...
            ,DataLabels{8},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})
        
        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,2),'r','LineWidth',2)
        legend({'Unfiltered', 'Comb MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{9},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})
    
        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{3},TeststoPlot(iTest)),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,Trig/10000,'b','LineWidth',2)
        hold on
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        
        plot(Time,vEMG(:,3),'r','LineWidth',2)
        legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})
        
        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,3),'r','LineWidth',2)
        legend({'Unfiltered', 'GS MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

        figure(10)
        subplot(2,1,1)
        plot(Time,Trig/10000,'b','LineWidth',2)
        hold on
        plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        
        plot(Time,vEMG_filtdropped(:,3),'r','LineWidth',2)
        legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})
        
        subplot(2,1,2)
        plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave_filtdropped(:,3),'r','LineWidth',2)
        legend({'Unfiltered', 'GS MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

    end
end


%% Plotting the dropped frames 
close all
iTrial=16;
TrialLabel=sprintf('Trial_%d',iTrial);
ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
FiltLabel="GS";
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    if ~S.(TestStruct).(ExpLabel).dropped 
                    figure(iTrial)
            subplot(length(TestFolders),1,iTest)
            text(0.3,0.5,'No Dropped Frames for This Trial')
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
        continue;
    else
       DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;   
        if isempty(DroppedFrames)
            figure(iTrial)
            subplot(length(TestFolders),1,iTest)
            text(0.3,0.5,'No Dropped Frames for This Trial')
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            continue;
        end

        DataLabels=S.(AnaLabel).AnaPar.DataLabels;
        DroppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG(1:end-3,:);
        FiltDropped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(:,DroppedFrames);
        clear lgd

        for iFrame=1:length(DroppedFrames)
            figure(iTrial)
            subplot(length(TestFolders),2,iTest*2-1)
            plot(DroppedEMG(:,iFrame),'LineWidth',2)
            hold on
            lgd(iFrame)=sprintf("Frame %d",DroppedFrames(iFrame));
            legend(lgd,'Location','NorthWest');
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            xlabel('Samples')
            ylabel('BP Filtered EMG')
            subplot(length(TestFolders),2,iTest*2)
            plot(FiltDropped(:,iFrame),'LineWidth',2)
            hold on 
            ttl=sprintf('Filtered Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            ylabel('BP Filtered and m-Wave Filtered')


            % ylim([3 -3]*10^-3)

        end
    end
end

%% Plotting EMG Features 
clc
AnaLabel=sprintf('%s_ana',TestFolders(1));
ExpLabel=S.(AnaLabel).AnaPar.ExpTable(1,:).('Occ');
TimeRange=[0.1 13];
TrialNum=[  4];
cm=lines(3);
FiltLabel="GS";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=ceil([stim_freq*TimeRange(1): stim_freq*TimeRange(2)]);

    for iTrial=1:length(TrialNum)
    
        TrialLabel=sprintf('Trial_%d',TrialNum(iTrial));
        AmpMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Amp_MAV_vEMG');
        AmpClipped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Clip_MAV_vEMG');   
        
        figure(1)
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,AmpMAV,'LineWidth',2,'DisplayName',sprintf('AmpMAV, Trial: %d',TrialNum(iTrial)),'Color',cm(1,:))
        hold on 
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,AmpClipped,'DisplayName',sprintf('AmpClipped, Trial: %d',TrialNum(iTrial)),'Color',cm(2,:))
        title(TestFolders(iTest))
    end
    
%     title(ttl)
    legend('Location','NorthWest')
    grid on
end
%% Plotting No Stim Trials

cm = lines(length(TestFolders));
sMVC=0;
FiltLabel="Unfilt";
Extrapolate=100;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');

    sMVCzero_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC');
    MAVMean_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean');
    vMVCVal_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC');
    
    % sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
    % MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps');
    % vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
    % Ind=(sMVCzero==sMVC);
    
    Amp_Modul_Mean=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean');
    
    Ind_Reps=(sMVCzero_Reps==sMVC);
    p1=polyfit(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),1);
    
    p2=polyfit(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),1);
    
    figure(1)
    subplot(2,1,1)
    plot(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p1,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    xlabel('% Voli. MVC')
    ylabel('Mean MAV')
    % xlim([0 50])
    
    subplot(2,1,2)
    plot(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p2,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    ylabel('% Amp modul')
    xlabel('% Voli. MVC')
    % xlim([0 50])
end



%% Comparing Theoretical and Actual MVC MAV 

NumofUpdate=4; % table mvc_table is updated 4 times inside loop below
NumofVariables=4;
mvc_table=strings(NumofUpdate*length(TestFolders),NumofVariables);
MVC_Percent=100;
sMVC_Ref=0;
AvgTimeMVC=2;
FiltLabel="Unfilt";
cm = lines(length(TestFolders));

for iTest=1:length(TestFolders)
    
    % MAV_MAX = avg MAV at the last 2 seconds
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    stim_freq=S.(TestLabel).ExpPar.stim_freq;

    AvgInd=stim_freq*AvgTimeMVC-10;

    TrialsMAV=zeros(1,S.(TestLabel).(ExpLabel).NumofTrials);
    TrialsAmp_MAV=zeros(1,S.(TestLabel).(ExpLabel).NumofTrials);
    for iTrial=1:S.(TestLabel).(ExpLabel).NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial );

        TrialsMAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
            .Feats(end-AvgInd:end,:).('MAV_vEMG'));
        
        TrialsAmp_MAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
            .AmpModulFeats(end-AvgInd:end,:).('Amp_MAV_vEMG'));
    end
    
    MAV_max=max(TrialsMAV);
    
    AmpModul_MAV_max=max(TrialsAmp_MAV);

    % Theoretical MAV_MAX
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    sMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC'));
    vMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC'));
    
    MAVMean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean'));
    Amp_Modul_Mean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean'));
    
    Ind=(sMVCLevs==sMVC_Ref);
    
    poly1=polyfit(vMVCLevs(Ind),MAVMean(Ind),1);
    poly2=polyfit(vMVCLevs(Ind),Amp_Modul_Mean(Ind),1);
    
    MAV_max_theo=polyval(poly1,MVC_Percent);
    AmpModul_max_theo=polyval(poly2,MVC_Percent);
%%-----> Plotting the MVC_ MAV levels

    figure(1) 
    subplot(2,1,1)
    semilogy(vMVCLevs(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    hold on 
    semilogy([unique(vMVCLevs(Ind)); 100],polyval(poly1,[unique(vMVCLevs(Ind)); 100]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest));
    semilogy(100, MAV_max, '*','Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    xlabel('% Voli. MVC')
    ylabel('Mean MAV')
    % xlim([0 50])
    
end
