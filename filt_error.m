%% Filtering Performance evaluation 
% # Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);
% # Mean MAV of RC Trials 
%%
clear MAV_Mean
exp_lbl='RC';
MeanTime=[8 10]; % Calculating the means at time [8 10]

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(exp_lbl);
    stim_freq=S.(TestLabel).ExpPar.stim_freq;

    MeanFrame=[MeanTime(1)*stim_freq MeanTime(2)*stim_freq];
    MeanRangeInd=[MeanFrame(1): MeanFrame(2)];
    FiltLabel="Unfilt";
    
    NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials';
    PercentMVCVals=S.(TestLabel).(ExpLabel).PercentMVC';
    PWPoints=S.(TestLabel).(ExpLabel).PWPoints';
    PWofTrials=S.(TestLabel).(ExpLabel).PWTrials';
    
    for iPW=1:length(PWPoints)
        
        Ind_PW(:,iPW)=PWofTrials==PWPoints(iPW);
        IndTrials(:,iPW)=find(Ind_PW(:,iPW)==1);
        
        for iTrial=1:length(IndTrials(:,iPW))
            TrialLabel=sprintf('Trial_%d',IndTrials(iTrial,iPW));
            MAV_Vals(iTrial,iPW)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel)...
                .(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
        end
    end
    
    MAV_Mean(:,iTest)=mean(MAV_Vals)';
    S.(AnaLabel).(ExpLabel).MAV_Vals=table(PWPoints, MAV_Vals(1,:)',...
        MAV_Vals(2,:)',MAV_Vals(3,:)',mean(MAV_Vals)',std(MAV_Vals)',...
        'VariableNames',["PW" "Rep_1" "Rep_2" "Rep_3" "Mean" "Std"]);
    
    S.(AnaLabel).(ExpLabel).MAV_Mean=array2table(MAV_Mean(:,iTest),...
        'VariableName',TestFolders(iTest));
    
end

% Update this to more general
RC_MAV=array2table(MAV_Mean,'VariableNames',TestFolders);


%Mean MAV of Occ Trials 
clc
exp_lbl='Occ';
VoliLevels=[10 20 30 40];
StimLevels=[ 0 10 20 30];
MeanTime=[ 11 14];
% Calculating the means at time [10 15]
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(exp_lbl);

    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    stim_freq=S.(TestLabel).ExpPar.FreqList(1);

    MeanFrame=[MeanTime(1)*stim_freq MeanTime(2)*stim_freq];
    MeanRangeInd=[MeanFrame(1):MeanFrame(2)];
    FiltLabel="Unfilt";
    c=1;
    clear MAV_Mean MAV_Mean_Reps TrialsInd
    for iVoli=1:length(VoliLevels)
        VoliLevel=VoliLevels(iVoli);

        for iStim=1:length(StimLevels)
            StimLevel=StimLevels(iStim);

            TrialsInd(:,c)= find_trialnum(VoliLevel, StimLevel, RepTableMat);
            
            for iTrial=1:length(TrialsInd(:,c))
                TrialLabel=sprintf('Trial_%d',TrialsInd(iTrial,c));

                MAV_Mean(TrialsInd(iTrial,c),1)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(MeanRangeInd).('MAV_vEMG'));
                MAV_Mean(TrialsInd(iTrial,c),2)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(MeanRangeInd).('MAV_vEMG'));
                MAV_Mean(TrialsInd(iTrial,c),3:5)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);
            end
            MAV_Mean_Reps(c,:) = [ c mean(MAV_Mean(TrialsInd(:,c),1))  std(MAV_Mean(TrialsInd(:,c),1)) RepTableMat(TrialsInd(iTrial,c),[4 5 ])];
            c=c+1;
        end
    end

    MAV_Mean=array2table(MAV_Mean,'VariableNames',["Mean", "Std" "Theo_Force","vMVC","sMVC"]);
    MAV_Mean_Reps=array2table(MAV_Mean_Reps,'VariableNames',["Unique_Trial" "Mean" "Std" "vMVC" "sMVC"]);
    
    S.(AnaLabel).(ExpLabel).MAV_Mean=MAV_Mean;
    S.(AnaLabel).(ExpLabel).MAV_Mean_Reps=MAV_Mean_Reps;
    S.(AnaLabel).(ExpLabel).RepTrialsInd=[array2table([TrialsInd'],...
        'VariableNames',["Rep1" "Rep2" "Rep3"]) MAV_Mean_Reps(:,'sMVC') MAV_Mean_Reps(:,'vMVC')];
end
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
    
    Ind=(sMVCzero==sMVC);
    Ind_Reps=(sMVCzero_Reps==sMVC);
    p=polyfit(vMVCVal(Ind),MAVMean(Ind),1);
    
    figure(10)
    plot(vMVCVal(Ind),MAVMean(Ind),'*','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot(vMVCVal(Ind),polyval(p,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
end

xlabel('% Voli. MVC')
ylabel('Avg MAV')
title(sprintf('Avg MAV (constant region) with Stim=%d%% MVC',sMVC))
grid
xlim([0 50])

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



