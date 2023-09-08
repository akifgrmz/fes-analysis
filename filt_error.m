%% Filtering Performance evaluation 
% # Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "apr20" "may19"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% Plotting

cm = lines(length(TestFolders));
sMVC=0;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
    
    sMVCzero_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
    MAVMean_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('Mean');
    vMVCVal_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('vMVC');
    
    sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC'); 
    MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean.('Mean');
    vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');
    
    Amp_Modul_Mean=S.(AnaLabel).(ExpLabel).Amp_Modul_Mean.('Mean');
    
    Ind=(sMVCzero==sMVC);
    Ind_Reps=(sMVCzero_Reps==sMVC);
    p1=polyfit(vMVCVal(Ind),MAVMean(Ind),1);

    p2=polyfit(vMVCVal(Ind),Amp_Modul_Mean(Ind),1);

    
    figure(1)
    subplot(2,1,1)
    plot(vMVCVal(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot(vMVCVal(Ind),polyval(p1,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    xlabel('% Voli. MVC')
    ylabel('Mean MAV')
    xlim([0 50])


    subplot(2,1,2)
    plot(vMVCVal(Ind),Amp_Modul_Mean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot(vMVCVal(Ind),polyval(p2,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    ylabel('% Amp modul')
    xlabel('% Voli. MVC')
    xlim([0 50])


end

% ylabel('Avg MAV')
% title(sprintf('Avg MAV (constant region) with Stim=%d%% MVC',sMVC))
% grid
%% Plotting Adaptive Modul
TestFolders="apr20";
AnaLabel=sprintf("%s_ana",TestFolders);

S = load_test(TestFolders,AnaLabel);
%%
PlotRange=[ 5 15];
vMVC=10;
sMVC=0;
FiltLabel="GS";
TestFolders="apr20";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
    
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=ceil([stim_freq*PlotRange(1): stim_freq*PlotRange(2)]);

%     ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);

    for iTrial=1:length(IndTrials)
        TrialLabel=sprintf("Trial_%d",IndTrials(iTrial));
        Adap_vEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Amp_MAV_vEMG');
        Clip_vEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Clip_MAV_vEMG');
        TrialNum=IndTrials(iTrial);
        cm=lines(length(IndTrials));

        MAV_vEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameInd,:).('MAV_vEMG');

        figure(1)    
        subplot(3,1,1)
        plot(FrameInd,MAV_vEMG,'DisplayName',sprintf('AmpMAV, Trial: %d',TrialNum),'Color',cm(iTrial,:))
        hold on
        ylabel('MAV')
        xlabel('Frames')
        title(sprintf('Filt: %s, Trial: %d, Stim: %d%% Voli: %d%%',...
        FiltLabel,TrialNum,RepTableMat(TrialNum,5),RepTableMat(TrialNum,4)));
        xlim([FrameInd(1) FrameInd(end)])
        grid on
        subplot(3,1,2)
        plot(FrameInd,Clip_vEMG,'DisplayName',sprintf('AmpClipped, Trial: %d',TrialNum),'Color',cm(iTrial,:))
        hold on
        ylabel('Norm. MAV')
        xlabel('Frames')
        title('Scaled MAV');
        xlim([FrameInd(1) FrameInd(end)])
        grid on


        subplot(3,1,3)
        plot(FrameInd,Adap_vEMG,'DisplayName',sprintf('AmpMAV, Trial: %d',TrialNum),'Color',cm(iTrial,:))
        hold on
        plot([175 350 525],[0 mean(Adap_vEMG(end-150:end)) mean(Adap_vEMG(end-150:end))],'Color',cm(iTrial,:))
        ylabel('Norm. MAV')
        xlabel('Frames')
        title('Amplitude Modulated MAV');
%         ylim([0  mean(Adap_vEMG(end-150:end))*1.5])
        xlim([FrameInd(1) FrameInd(end)])
        grid on

        figure(2)    
        

        

    end
end



%%


figure
for iTest=1:length(TestFolders)
    plot(linspace(10,30,7),RC_MAV.(TestFolders(iTest)),'DisplayName',TestFolders(iTest))
    hold on
end

xlabel('Stim MVC (%)')
ylabel('MWave MAV')
title('Average MAV from RC Trials')
grid on 
legend()

%%
% # Comparing MWave and Voli MAV

for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
    sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
    OccMeans=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps(sMVC==0,:).('Mean');    
    OccStds=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps(sMVC==0,:).('Std');
    MVCRange=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps(sMVC==0,:).('vMVC');
    
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('RC');
    RCMeans=S.(AnaLabel).(ExpLabel).MAV_Vals.('Mean');
    RCStds=S.(AnaLabel).(ExpLabel).MAV_Vals.('Std');
    PWPoints=S.(AnaLabel).(ExpLabel).MAV_Vals.('PW');

    figure(1)
    subplot(length(TestFolders),2,2*iTest-1)
    bar(MVCRange,OccMeans)
    hold  
    er = errorbar(MVCRange,OccMeans,-OccStds,+OccStds);    
    er.Color = [0 0 0];
    er.LineStyle = 'none'; 
    ylabel(TestFolders{iTest},'fontweight','bold','fontsize',14)
    ylim([ 0 4*10^-4])
    xlabel("Percent MVC")
    grid

    subplot(length(TestFolders),2,2*iTest)
    bar(PWPoints,RCMeans)
    hold
    er = errorbar(PWPoints,RCMeans,-RCStds,+RCStds);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylim([ 0 4*10^-4])
    xlabel("PW")
    grid
end

subplot(length(TestFolders),2,1)
title("vEMG MAV (No Stimulation)",'fontweight','bold','fontsize',14)

subplot(length(TestFolders),2,2)
title("M-Wave MAV (No Voli. Effort )",'fontweight','bold','fontsize',14)

%% 
% # Evaluate the filter accuracy
% NRMSE_s= sqrt(sum((MAV_voli-MAV_filt)^2)/T)/std(MAV_voli)
% NRMSE_v= sqrt(sum((MAV_voli_1-MAV_voli_2)^2)/T)/std(MAV_voli_1)

c=1;

TimeRange=[5 10] ;
ExpstoAna= ["Occ"];
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" ];
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    clear ExpLabels
    for iExp=1:length(ExpstoAna)
        ExpLabels(iExp)=S.(AnaLabel).AnaPar.ExpTable.(ExpstoAna(iExp));
    end

    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameRange=[TimeRange(1)*stim_freq TimeRange(2)*stim_freq];
    FrameRangeInd=[FrameRange(1): FrameRange(2)];
    
    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels(iExp);
        
        NumofTrial=S.(TestLabel).(ExpLabel).NumofTrials;
        
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat,...
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);

        ref_count=1;
        for iTrial=1:NumofTrial
            TrialLabel=sprintf("Trial_%d",iTrial);
            FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
            
            BegofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);
            
            MVC_Voli=RepMatTable(iTrial,:).('MVC_Voli');
            MVC_Stim=RepMatTable(iTrial,:).('MVC_Stim');
            
            RefIndTrial=find_trialnum(MVC_Voli,0,S.(TestLabel).(ExpLabel).RepTableMat);

            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels{iFilt};

                MAV_filt=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameRangeInd,:).('MAV_vEMG');
                clear MAV_ref
                for iRefTrial=1:length(RefIndTrial)
                    TrialLabel2=sprintf("Trial_%d",RefIndTrial(iRefTrial));

                    MAV_ref(:,iRefTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel2).Unfilt.Feats(FrameRangeInd,:).('MAV_vEMG');
                    rmse(iRefTrial,iFilt)=sqrt( sum((MAV_filt-MAV_ref(:,iRefTrial)).^2)/length(MAV_filt));%/mean(MAV_ref(:,iRefTrial));
                end
            end
            RMSE=table(rmse(:,1),rmse(:,2),rmse(:,3),RefIndTrial,[MVC_Voli;MVC_Voli;MVC_Voli], ...
                [MVC_Stim;MVC_Stim;MVC_Stim],'VariableNames',[ FiltLabels(1) FiltLabels(2) FiltLabels(3)...
                "Ref_Trial" "MVC_Voli" "MVC_Stim"]);
            
            S.(AnaLabel).(ExpLabel).(TrialLabel).RMSE=RMSE;
            clear RMSE
            
        end
    end
end

%% Plotting NRMSE
clc
VoliLevels=[10 20 30 40];
StimLevels=[ 0 10 20 30];

TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" ];
AnaLabel=sprintf("%s_ana",TestFolders{1});
FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');

for iVoli=1:length(VoliLevels)
    VoliLevel=VoliLevels(iVoli);

    for iStim=1:length(StimLevels)
        StimLevel=StimLevels(iStim);
        
        NRMSE=zeros(0,length(FiltLabels));
        
        for iTest=1:length(TestFolders)
            TestLabel=sprintf("%s_test",TestFolders{iTest});
            AnaLabel=sprintf("%s_ana",TestFolders{iTest});

            RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
            stim_freq=S.(TestLabel).ExpPar.stim_freq;

            TrialsInd= find_trialnum(VoliLevel, StimLevel, RepTableMat);
            for iTrial=1:length(TrialsInd)
                TrialLabel=sprintf("Trial_%d",TrialsInd(iTrial));
                
                NRMSE=[NRMSE; 
                    [S.(AnaLabel).(ExpLabel).(TrialLabel).RMSE.(FiltLabels(1))'
                     S.(AnaLabel).(ExpLabel).(TrialLabel).RMSE.(FiltLabels(2))'
                     S.(AnaLabel).(ExpLabel).(TrialLabel).RMSE.(FiltLabels(3))']'  ];
            end
        end
        if StimLevel==0
            NRMSE(:,3)=[];
            NRMSE(:,2)=[];
            FiltLabelsPlot=FiltLabels(1);
            NRMSE=NRMSE(NRMSE~=0);
        else
            FiltLabelsPlot=FiltLabels ;

        end
        
        figure(1)
        subplot(length(StimLevels),length(VoliLevels),(iVoli-1)*(length(VoliLevels))+iStim)
        bar(categorical(FiltLabelsPlot),mean(NRMSE) )
        [r,c]=size(NRMSE);
        hold 
        er = errorbar(categorical(FiltLabelsPlot),mean(NRMSE),-std(NRMSE),+std(NRMSE));    
        er.Color = [0 0 0];
        er.LineStyle = 'none'; 
        ylim([ 0 2*10^-4])
        grid
        subplot(length(StimLevels),length(VoliLevels),iStim)
        title(sprintf("Stim: %d%%, n= %d(each bar)",StimLevel,r),'fontweight','bold','fontsize',14)

    end
    subplot(length(StimLevels),length(VoliLevels),(iVoli-1)*length(StimLevels)+1)
    ylabel(sprintf('Voli: %d%%',VoliLevel),'fontweight','bold','fontsize',14)

end



