%% Plotting for EMBC 
clc
clear all
TestFolders=[ "jan7"];

for iTest=1:length(TestFolders)
    TestFiles{iTest}=sprintf("%s_ana",TestFolders{iTest});
end

AnaStruct=sprintf("%s_ana",TestFolders);
S = load_test(TestFolders,TestFiles);
%% 
close all
TestStruct=sprintf("%s_test",TestFolders);
AmpGain=S.(TestStruct).ExpPar.AmpGain;
Lbl='Occ';
ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
PlotTrial=[15 19];
PlotFrame=floor([ 486 490]);
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
    plot(Time,reshape(MWave,1,[]),'k','LineWidth',2)
    hold on
    plot(Time,reshape(vEMG,1,[]),'r','LineWidth',2)
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
PlotTrial=[ 7 7];
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
    PWofTrials=S.(TestStruct).(ExpLabel).PWTrials;
    PWPoints=S.(TestStruct).(ExpLabel).PWPoints;
    PWProfile=S.(TestStruct).(ExpLabel).PWProfile;
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
    
    plot(Target(1,1:4),Target(2,1:4)*RepTableMat(IndTrials(iTrial),1),'b','LineWidth',2)
    plot(T,PW/30,'r','LineWidth',2)

    lgd{end+1}="Target Line";
    lgd{end+1}="PW (a.u.)";

    legend(lgd)
    grid on
    title('Force Measurement Superimposed to Visual Target Line')
    xlabel('Time(s)')
    ylabel('Force(N)')
    
end
%% statistical tests results 


%% Plot Dropped vs Filtered
clear IndTrials vMVC sMVC IndS IndV DroppedFramesIndexFixed DroppedFrames UnfiltFeat FiltMAV vEMGFeat DroppedFeat DroppedFramesIndexFixed DroppedFrameInd x
close all 
clc

MarginFromDropped=5;  % frames
PlotRange1=[ 5.4 10]; 
PlotRange2=[ 10 15.2]; 
PlotVoli=1; VoliMVCLevels=[10 20 30 40 ];
PlotStim=2; StimMVCLevels=[ 0 10 20 30];
PlotRangeFrames1=PlotRange1*stim_freq;  
PlotRangeFrames2=PlotRange2*stim_freq;        
PlotRangeInd=[PlotRangeFrames1(1):PlotRangeFrames1(2),PlotRangeFrames2(1):PlotRangeFrames2(2)];
PreDroppedRange= [10.8 11 ];
PostDroppedRange=[11.2 11.5]; %% based on the dropped indices 
DroppedRange1=[5.1 10.8]*stim_freq;
DroppedRange2=[16.1 18.9]*stim_freq;
boolean DroppedFrameInd;
TrialColor={'r', 'b', 'k'};
iFeat=1;
FeatInd=[ 1 1];
exp_lbl='Occ';


for iTest=1:1
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
    FeatLabel=S.(AnaStruct).AnaPar.FeatLabels{FeatInd(iFeat)};
    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
    sMVC=StimMVCLevels(PlotStim);
    vMVC=VoliMVCLevels(PlotVoli);
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

    clear FiltRelatedFramesMean FiltRelatedFramesSD UnfiltRelatedFramesMean UnfiltRelatedFramesSD  PercentDrops
    clear UnfiltRelatedFrames FiltRelatedFrames
    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;
    TargetFrames=mean(S.(AnaStruct).(ExpLabel).TargetFrames);
    for iFilt=2:length(S.(AnaStruct).AnaPar.FiltLabels)
        FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};
        vEMGLabel=sprintf('Norm_%s_vEMG',FeatLabel);
        DroppedFeatLabel=sprintf('Norm_%s_vEMG',FeatLabel);

        for iPlot=1:length(IndTrials)
            TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));

            DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
            DroppedFrameInd1=PlotRangeFrames1(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames1(2);
            DroppedFrameInd2=PlotRangeFrames2(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames2(2);
            DroppedFrameInd=DroppedFrameInd1+DroppedFrameInd2;
            DroppedFrameInd=logical(DroppedFrameInd);
            
            vEMGFeat(iPlot,:)=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel)(PlotRangeInd);
            
%             vEMGFeat(iPlot,:)=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,PlotRangeInd);
            DroppedFramesFeat(iPlot,:)=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel)(DroppedFrameInd);
            
            UnfiltFeat(iPlot,:)=S.(AnaStruct).(ExpLabel).(TrialLabel).Unfilt.Feats.(vEMGLabel)(PlotRangeInd);
%             FiltMAV(iPlot,:)=P.(ExpLabel).(TrialLabel).(FiltLabel).FiltMAV(1,:); % 1 vEMG 2 MWaves
            font=14;
            figure(1)
            set(gca,'FontSize',font)
            subplot(2,1,iFilt-1)
            for iDropped=1:length(DroppedFrames(DroppedFrameInd))
                DroppedLabel=sprintf('Dropped_%d',iDropped);
                
                DroppedInd=DroppedFrames(iDropped);
                FiltRelatedFramesMean(iPlot,iDropped)=mean(S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                FiltRelatedFramesSD(iPlot,iDropped)=std(S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                UnfiltRelatedFramesMean(iPlot,iDropped)=mean(S.(AnaStruct).(ExpLabel).(TrialLabel).Unfilt.Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                UnfiltRelatedFramesSD(iPlot,iDropped)=std(S.(AnaStruct).(ExpLabel).(TrialLabel).Unfilt.Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                UnfiltRelatedFrames.(DroppedLabel){iPlot}=S.(AnaStruct).(ExpLabel).(TrialLabel).Unfilt.Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);
                FiltRelatedFrames.(DroppedLabel){iPlot}=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);
                g1=S.(AnaStruct).(ExpLabel).(TrialLabel).Unfilt.Feats.(vEMGLabel)...
                    (-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);
                g2=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel)(DroppedFrameInd);
                PlotInd=-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped;
                OnesInd=ones(length(PlotInd),1);
                plot(PlotInd,OnesInd*UnfiltRelatedFramesMean(iPlot,iDropped),...
                    'LineWidth',3,'Color','k')
                hold on
                plot(PlotInd,OnesInd*FiltRelatedFramesMean(iPlot,iDropped),...
                    'LineWidth',3,'Color','r')
            end

            plot(PlotRangeInd,vEMGFeat(iPlot,:),'.', 'Color','r','LineWidth',0.05)
            plot(PlotRangeInd,UnfiltFeat(iPlot,:),'.','Color','k','LineWidth',0.05)
            plot(DroppedFrames(DroppedFrameInd),DroppedFramesFeat(iPlot,:),'o','Color','b','LineWidth',1.5)
            TrialString=num2str(IndTrials'); 
            ttl = sprintf('P-values(95%% CL) with %s filtered MAV, Stim: %d%%, Voli: %d%%, Trials: %s',FiltLabel,sMVC,vMVC,TrialString);
            title(ttl);
            str=sprintf("%s Filtered Mean Around Dropped",FiltLabel);
            legend(["Unfiltered Mean Around Dropped"; str] ,'Location','northwest')
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'k', 'DisplayName', sprintf('Unfiltered MAV'))
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'r', 'DisplayName', sprintf('%s Filtered MAV',FiltLabel))
            plot([NaN NaN], [NaN NaN],'o', 'Color', 'b', 'DisplayName', sprintf('Dropped Frames'))
            plot(PlotRangeInd,(4+4*TargetFrames(PlotRangeInd))/100000,'Color', 'r','DisplayName',sprintf('Norm. Target Line'));
            ylabel('EMG MAV');
            xlabel('Frames');
        end
        
        PercentDrops=(UnfiltRelatedFramesMean-DroppedFramesFeat)./UnfiltRelatedFramesMean*100;
        figure
        gLabels={'Unfilt' 'Filtered' 'Dropped'};
        for iDropped=1:length(DroppedFrames(DroppedFrameInd))
            for iPlot=1:length(IndTrials)
                DroppedLabel=sprintf('Dropped_%d',iDropped);

                Unfilts(iPlot,:)=UnfiltRelatedFrames.(DroppedLabel){iPlot};
                y1=reshape(Unfilts,1,[]);
                g1=cell(1,length(y1));
                g1(:)=gLabels(1);
                Filts(iPlot,:)=FiltRelatedFrames.(DroppedLabel){iPlot};
                y2=reshape(Filts,1,[]);
                g2=cell(1,length(y2));
                g2(:)=gLabels(2);
                y3(iPlot)=DroppedFramesFeat(iPlot,iDropped);
                g3=cell(1,length(y3));
                g3(:)=gLabels(3);
                g=[ g1 g2 g3];
                y=[y1 y2 y3];

            end
            
            [p(iDropped),t,stats]= anova1(y,g,'off');
            subplot(1,length(DroppedFrames(DroppedFrameInd)),iDropped)
            results{iDropped}=multcompare(stats);
            title([])
            ylabel([])
            xlabel([])

        end
                    figure(1)
            subplot(2,1,iFilt-1)
        for iDropped=1:length(DroppedFrames(DroppedFrameInd))
            DroppedInd=DroppedFrames(iDropped);

            p_vals=results{iDropped};
            text(DroppedInd-10,2.2*max(UnfiltRelatedFramesMean(1,:)),sprintf('p: %.3f',p_vals(1,6)), 'FontSize', font)
            text(DroppedInd-10,2*max(UnfiltRelatedFramesMean(1,:)),sprintf('p: %.3f',p_vals(2,6)), 'FontSize', font)
            text(DroppedInd-10,1.8*max(UnfiltRelatedFramesMean(1,:)),sprintf('p: %.3f',p_vals(3,6)), 'FontSize', font)
        end
        
        text(155,2.2*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(1,1)},gLabels{p_vals(1,2)}), 'FontSize', font)
        text(155,2.0*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(2,1)},gLabels{p_vals(2,2)}), 'FontSize', font)
        text(155,1.8*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(3,1)},gLabels{p_vals(3,2)}), 'FontSize', font)
        ylim([0 2.6*max(UnfiltRelatedFramesMean(1,:))])
                         

        for iPlot=1:length(IndTrials)
            TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));
            
            MeanDrop=mean( PercentDrops);
            StdDrop=std( PercentDrops);
        end
    end
end




%%















