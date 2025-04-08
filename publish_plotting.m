%% Load 
TestFolders=["feb28_24" "feb29_24" "mar20_24" ]; % these tests need manu_additoins
textt="%s_ana-nov17-";
TestFiles=compose(textt,[TestFolders']);
S = load_test(TestFolders,TestFiles);
%% EMG plotting
close all
Lbl='Occ';
PlotTrial=[ 10 11];
PlotTime=[ 12 13];
% PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);
PlotFrame=floor([ 278 282]);

TeststoPlot=TestFolders(end);

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

        % figure(10)
        % subplot(2,1,1)
        % plot(Time,Trig/10000,'b','LineWidth',2)
        % hold on
        % plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        % 
        % plot(Time,vEMG_filtdropped(:,3),'r','LineWidth',2)
        % legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
        % ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
        %     DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        % title(ttl);
        % xlabel(DataLabels{5})
        % ylabel(DataLabels{1})
        % 
        % subplot(2,1,2)
        % plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        % hold
        % plot(Time,MWave_filtdropped(:,3),'r','LineWidth',2)
        % legend({'Unfiltered', 'GS MWaves'})
        % ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
        %     DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        % title(ttl);
        % xlabel(DataLabels{5})
        % ylabel(DataLabels{1})

    end
end

%% Effort, Dropped, EMG

close all
c=[300 300 450];
clc
lbl='Occ';
PlotVoli=3;
PlotStim=3;
cm=lines(length(TestFolders));
TimeRange=[1 17];
FiltLabel="GS";
DroppedFiltLabel="Unfilt";

for iTest=1:1
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    % TestName=TestFolders(iTest);
    TestName=S.(AnaLabel).AnaPar.TestName;
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    MVC=S.(TestLabel).(ExpLabel).MVC;

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(lbl);
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
    sMVC=StimMVCLevels(PlotStim);
    vMVC=VoliMVCLevels(PlotVoli);
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

    StimOff=S.(TestLabel).(ExpLabel).StimProfile;
    PlotRange=[ StimOff-10 StimOff+1];
    PreStimOffRange= [StimOff-1 StimOff ];
    PostStimOffRange=[StimOff+.2 StimOff+.3];
    ttl=sprintf('BandPass Filtered MAV Compared to Dropped Frames (Stim: %d%%, Voli: %d%%)',sMVC,vMVC);
    EffortLabel=S.(TestLabel).ExpPar.EffortType;
    for iTrial=1:length(IndTrials)
    IndTrials(iTrial)
        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        EffortMea=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortLabel);
        MAV=S.(AnaLabel).(ExpLabel).(TrialLabel).Unfilt.Feats.('MAV_vEMG');

        FiltMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats.('FiltMAV_vEMG');
        AmpModulMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Amp_MAV_vEMG');

        MAVInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd;
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time'); 
        TimeInd=Time>TimeRange(1) & Time<TimeRange(2);

        DroppedMAV=[];
        DroppedInd=[];
        DroppedFeatFilt=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat.('Filt');
        if sMVC~=0
            DroppedMAV=str2double(S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(DroppedFeatFilt==DroppedFiltLabel,:).('MAV'));
            DroppedInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
        end
        
        % f=figure(1);
        % f.Position = [100 300 2000 800];    
        % subplot(2,1,1)
        % plot(Time(TimeInd),EffortMea(TimeInd),'Color',cm(iTest,:),'DisplayName',TestName,'LineWidth',2)
        % hold on
        % subplot(2,1,2)
        % plot(Time(TimeInd),EffortMea(TimeInd)/MVC*100,'Color',cm(iTest,:),'DisplayName',TestName,'LineWidth',2)
        % hold on

        f=figure(2);
        f.Position = [100 300 2000 800];          
        subplot(2,1,1)
        plot(MAVInd,MAV,'o','Color',cm(iTest,:),'LineWidth',1,'DisplayName',sprintf('MAV Before Filtering'))%sprintf('MAV (%s)',TestName)
        hold on
        plot(DroppedInd,DroppedMAV,'*','Color',cm(iTest,:),'LineWidth',3,'DisplayName',sprintf('Dropped MAV'))%sprintf('Dropped MAV (%s)',TestName))
                ylim([0 2*10^-4])

        subplot(2,1,2)
        % plot(MAVInd,FiltMAV,'Color',cm(iTest,:),'DisplayName',TestName)
        plot(MAVInd,AmpModulMAV/c(iTest),'Color',cm(iTest,:),'LineWidth',2.5,'DisplayName',sprintf('GS filtered MAV')) %sprintf('GS filtered MAV (%s)',TestName))
        hold on
        plot(DroppedInd,DroppedMAV,'*','Color',cm(iTest,:),'LineWidth',3,'DisplayName',sprintf('Dropped MAV'))%sprintf('Dropped MAV (%s)',TestName))
        ylim([0 2*10^-4])
        % subplot(3,1,3)
        % plot(MAVInd,AmpModulMAV,'Color',cm(iTest,:),'DisplayName',TestName)
        % hold on
%         plot(DroppedInd,DroppedMAV/MVC*100,'*','Color',cm(iTest,:),'DisplayName',TestName)

    end
end

% figure(1)
% subplot(2,1,1)
% grid on
% legend()
% title(ttl)
% ylabel("Measured Effort")
% subplot(2,1,2)
% grid on
% legend()
% ylabel("Percent MVC")
% title(ttl)

figure(2)
subplot(2,1,1)
ylabel("MAV")
legend()
% subplot(3,1,2)
ylabel("BandPass filtered MAV")
legend()
title(ttl)
xlabel("Frames")
ylim([0 0.0002])

subplot(2,1,2)
title("GS Filtered MAV Compared to Dropped Frames")
ylabel("MAV")
ylim([0 0.00005])
xlabel("Frames")
legend()

% ylim([0 0.001])
% subplot(3,1,3)
% grid on
% ylabel("Norm. MAV")
% % ylim([0 .3])