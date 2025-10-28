%% Occlusion Paper Figures 
% # Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "apr20" "may19"];
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24" "aug29_24" "sep3_24" "sep4_24" "sep6_24"];
TestFolders=[ "jan7" "jan11" "jan12"  "aug26_24" "aug29_24" "sep3_24" "sep4_24" "sep6_24"];
TestFolders=["jan7" "jan11" "jan12"  "sep3_24" "sep4_24" "sep6_24" "oct17_24" "oct18_24"];


TestFolders=["jan7" "jan11" "jan12"  "sep3_24" "sep4_24" "sep6_24" "oct17_24" "oct18_24"];

textt="%s_ana-oct19-";
TestFiles=compose(textt,[TestFolders']);
S = load_test(TestFolders,TestFiles);

%% Plots after filtering and preproccess 
clc
close all
Trials=[ 1: 10]; 
TimeRange=[1 12];

% TestLabel=sprintf('%s_test',FolderNames);
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
        % TimeRange=S.(TestLabel).(ExpLabel).StimRange;
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
        TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);

        % RawEMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('EMG');
        BlankedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('BlankedEMG');

        Trigger=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('Trigger');
        EffortMeas=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).(EffortType);
        PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('PW');

        f=figure(iTest);
        f.Position = [100 100 1400 800];                
        subplot(2,length(Trials),iTrial)
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
%%

%% Plotting, Debugging after Filtering
Lbl='Occ';
PlotTrial=[ 12 12];
PlotTime=[ 12 13];
% PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);
PlotFrame=floor([ 261 265]);
% DroppedFrames=264*98:265*98;

TeststoPlot=TestFolders(end);

for iTest=1:length(TeststoPlot)
    TestLabel=sprintf("%s_test",TeststoPlot(iTest));
    AnaLabel=sprintf("%s_ana",TeststoPlot(iTest));
    DataLabels=S.(AnaLabel).AnaPar.DataLabels;
    TestName=S.(AnaLabel).AnaPar.TestName;
    ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).(Lbl));
    cm=lines(5);
    % PlotFrame=floor([ 300 305]); % first frame is skipped
    % ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    % ExpstoAna=ExpLabels([2,3,4,5]);
    FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
    
    clear EMG EMGTime EMGFrames MWaveFrames vEMG MWave MWaveFrames_filtdropped vEMGFrames_filtdropped vEMG_filtdropped MWave_filtdropped MWave_withdropped vEMGFrames_withdropped MWaveFrames_withdropped vEMG_withdropped
    % vEMGFrames_filtdropped=zeros(FrameLength,PlotFrame(2)-PlotFrame(1)+1,length(FiltLabels));
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameLength=S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength;
        DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
        

        DroppedPlot=DroppedFrames((PlotFrame(1)<=DroppedFrames) & DroppedFrames<=PlotFrame(2));
        DroppedFrameSamp=DroppedPlot*98:(DroppedPlot+1)*98;


        TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
        EMG=S.(AnaLabel).(ExpLabel).(TrialLabel).EMGFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));

        TriggerFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TriggerFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
    
        Time=reshape(TimeFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        Trig=reshape(TriggerFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        EMGTime(iTest,:)=reshape(EMG,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

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
    
            vEMGFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
             
            vEMG_withdropped(:,iFilt)=reshape(vEMGFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave_withdropped(:,iFilt)=reshape(MWaveFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            
            droppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped(1:FrameLength,DroppedPlot);
            droppedTime=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,DroppedPlot);
        end
    end
end

figure(100)
clf
set(gcf,'Position', [1800, 400, 600, 300]);
hold on


% plot(Time,EMGTime,'Color',cm(1,:),'LineWidth',2,'DisplayName',strcat('Before blanking'));
plot(Time,vEMG_withdropped(:,1),'Color',cm(1,:),'LineWidth',2,'DisplayName',strcat('Unfiltered EMG'));
plot(Time,MWave_withdropped(:,2),'Color',cm(2,:),'LineStyle','-','LineWidth',2,'DisplayName',strcat('Est. mWave'));
plot(Time,vEMG_withdropped(:,3),'Color',cm(3,:),'LineStyle','-','LineWidth',2,'DisplayName',strcat('GS Filtered vEMG'));
plot(droppedTime,droppedEMG,'Color','k','LineWidth',2,'DisplayName',strcat('Dropped Frame'));
xlim([Time(1) Time(end)])
xlabel('Time(s)')
ylabel('EMG (mV)')
grid on 
legend
grid on
legend('Location','SouthEast');
set(gca, 'FontSize', 14);  % Optional: makes axis tick labels larger too


% Add text with arrow pointing to a specific data point
arrow1beg = [7.29057 -0.000241404];
arrow1end = [7.3 -0.00031];

arrow2beg = [7.31171 0];
arrow2end = [7.31 0.00016];

arrow3beg = [7.33171 0.000027];
arrow3end = [7.33171 0.00021];

arrow1length=[0.1 -10^-5 ]; % normalized
x2_data = 7.3117;
y2_data = 0;
arrow2length=[-0.05 0.2]; % normalized

xlims = xlim;
ylims = ylim;
beg1norm = data2norm(arrow1beg,xlims,ylims);
end1norm = data2norm(arrow1end,xlims,ylims);

beg2norm = data2norm(arrow2beg,xlims,ylims);
end2norm = data2norm(arrow2end,xlims,ylims);


beg3norm = data2norm(arrow3beg,xlims,ylims);
end3norm = data2norm(arrow3end,xlims,ylims);

% [norm_x2, norm_y2] = data2norm(x2_data, y2_data,xlims,ylims);

% Add arrow from text to point
annotation('arrow', [beg1norm(1), end1norm(1)], [beg1norm(2), end1norm(2)], 'LineWidth', 1.5);

text(arrow1end(1)+0.01, arrow1end(2), 'mWave', ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

annotation('arrow', [beg2norm(1), end2norm(1)], [beg2norm(2), end2norm(2)], 'LineWidth', 1.5);

text(arrow2end(1), arrow2end(2)+0.00005, ["Blanked";"Period  "], ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');

annotation('arrow', [beg3norm(1), end3norm(1)], [beg3norm(2), end3norm(2)], 'LineWidth', 1.5);

text(arrow3end(1)+0.005, arrow3end(2)+0.00004, ["vEMG"], ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top');


function norm = data2norm(xy,xlims,ylims)
    % Get axes position in figure (normalized units)
    ax_pos = get(gca, 'Position');  % [left, bottom, width, height]
    
    % Get axes limits
    % xlims = xlim;
    % ylims = ylim;
    x_data=xy(1);
    y_data=xy(2);
    % Convert data coordinates to normalized axes coordinates (0 to 1 within axes)
    x_norm_ax = (x_data - xlims(1)) / (xlims(2) - xlims(1));
    y_norm_ax = (y_data - ylims(1)) / (ylims(2) - ylims(1));
    
    % Convert to normalized figure coordinates
    norm_x = ax_pos(1) + x_norm_ax * ax_pos(3);
    norm_y = ax_pos(2) + y_norm_ax * ax_pos(4);

    norm=[norm_x norm_y] ;
end


%% 
% plot(Time,MWave_withdropped(:,iFilt),'Color',cm(iFilt,:),'LineStyle','--','LineWidth',2,'DisplayName',strcat('mWaves',' (',TrialLabel,',',FiltLabel));

% figure('Name',sprintf('%s Results, Test: %s',FiltLabels{2},TeststoPlot(iTest)),'NumberTitle','off')
% subplot(2,1,1)
% plot(Time,Trig/10000,'b','LineWidth',2)
% hold on
% plot(Time,vEMG(:,1),'k','LineWidth',2)
% 
% plot(Time,vEMG(:,2),'r','LineWidth',2)
% legend({'Trigger(a.u)', 'Unfiltered', 'Comb vEMG'})
% ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d'...
%     ,DataLabels{8},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
% title(ttl);
% xlabel(DataLabels{5})
% ylabel(DataLabels{1})
% 
% subplot(2,1,2)
% plot(Time,vEMG(:,1),'k','LineWidth',2)
% hold
% plot(Time,MWave(:,2),'r','LineWidth',2)
% legend({'Unfiltered', 'Comb MWaves'})
% ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%     DataLabels{9},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
% title(ttl);
% xlabel(DataLabels{5})
% ylabel(DataLabels{1})
% 
% figure('Name',sprintf('%s Results, Test: %s',FiltLabels{3},TeststoPlot(iTest)),'NumberTitle','off')
% subplot(2,1,1)
% plot(Time,Trig/10000,'b','LineWidth',2)
% hold on
% plot(Time,vEMG(:,1),'k','LineWidth',2)
% 
% plot(Time,vEMG(:,3),'r','LineWidth',2)
% legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
% ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%     DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
% title(ttl);
% xlabel(DataLabels{5})
% ylabel(DataLabels{1})
% 
% subplot(2,1,2)
% plot(Time,vEMG(:,1),'k','LineWidth',2)
% hold
% plot(Time,MWave(:,3),'r','LineWidth',2)
% legend({'Unfiltered', 'GS MWaves'})
% ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%     DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
% title(ttl);
% xlabel(DataLabels{5})
% ylabel(DataLabels{1})
% 

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
