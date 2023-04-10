%% Filtering Performance evaluation 
% # Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16"];

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
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);
    stim_freq=S.(TestStruct).ExpPar.stim_freq;

    MeanFrame=[MeanTime(1)*stim_freq MeanTime(2)*stim_freq];
    MeanRangeInd=[MeanFrame(1): MeanFrame(2)];
    FiltLabel="Unfilt";
    
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials';
    PercentMVCVals=S.(TestStruct).(ExpLabel).PercentMVC';
    PWPoints=S.(TestStruct).(ExpLabel).PWPoints';
    PWofTrials=S.(TestStruct).(ExpLabel).PWTrials';
    
    for iPW=1:length(PWPoints)
        
        Ind_PW(:,iPW)=PWofTrials==PWPoints(iPW);
        IndTrials(:,iPW)=find(Ind_PW(:,iPW)==1);
        
        for iTrial=1:length(IndTrials(:,iPW))
            TrialLabel=sprintf('Trial_%d',IndTrials(iTrial,iPW));
            MAV_Vals(iTrial,iPW)=mean(S.(AnaStruct).(ExpLabel).(TrialLabel)...
                .(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
        end
    end
    
    MAV_Mean(:,iTest)=mean(MAV_Vals)';
    S.(AnaStruct).(ExpLabel).MAV_Vals=table(PWPoints, MAV_Vals(1,:)',MAV_Vals(2,:)',MAV_Vals(3,:)',mean(MAV_Vals)',std(MAV_Vals)',...
        'VariableNames',["PW" "Rep_1" "Rep_2" "Rep_3" "Mean" "Std"]);
    
    S.(AnaStruct).(ExpLabel).MAV_Mean=array2table(MAV_Mean(:,iTest),...
        'VariableName',TestFolders(iTest));
    
end

% Update this to more general
RC_MAV=array2table(MAV_Mean,'VariableNames',TestFolders);


%
%Mean MAV of Occ Trials 
clc
exp_lbl='Occ';
VoliLevels=[10 20 30 40];
StimLevels=[ 0 10 20 30];
MeanTime=[ 11 14];
% Calculating the means at time [10 15]
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);

    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
    stim_freq=S.(TestStruct).ExpPar.FreqList(1);

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

                MAV_Mean(TrialsInd(iTrial,c),1)=mean(S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
                MAV_Mean(TrialsInd(iTrial,c),2)=std(S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
                MAV_Mean(TrialsInd(iTrial,c),3:5)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);
            end
            MAV_Mean_Reps(c,:) = [ c mean(MAV_Mean(TrialsInd(:,c),1))  std(MAV_Mean(TrialsInd(:,c),1)) RepTableMat(TrialsInd(iTrial,c),[4 5 ])];
            c=c+1;
        end
    end

    MAV_Mean=array2table(MAV_Mean,'VariableNames',["Mean", "Std" "Theo_Force","vMVC","sMVC"]);
    MAV_Mean_Reps=array2table(MAV_Mean_Reps,'VariableNames',["Unique_Trial" "Mean" "Std" "vMVC" "sMVC"]);
    
    S.(AnaStruct).(ExpLabel).MAV_Mean=MAV_Mean;
    S.(AnaStruct).(ExpLabel).MAV_Mean_Reps=MAV_Mean_Reps;
    S.(AnaStruct).(ExpLabel).RepTrialsInd=[array2table([TrialsInd'],...
        'VariableNames',["Rep1" "Rep2" "Rep3"]) MAV_Mean_Reps(:,'sMVC') MAV_Mean_Reps(:,'vMVC')];
end
%% Plotting
cm = lines(length(TestFolders));
ExpLabel=S.(AnaStruct).AnaPar.ExpTable.('Occ');

sMVC=0;
for iTest=1:length(TestFolders)
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    sMVCzero_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('sMVC');
    MAVMean_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('Mean');
    vMVCVal_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('vMVC');
    
    sMVCzero=S.(AnaStruct).(ExpLabel).MAV_Mean.('sMVC'); 
    MAVMean=S.(AnaStruct).(ExpLabel).MAV_Mean.('Mean');
    vMVCVal=S.(AnaStruct).(ExpLabel).MAV_Mean.('vMVC');
    
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
title("vEMG MAV (No Stimulation Trials)",'fontweight','bold','fontsize',14)

subplot(length(TestFolders),2,2)

title("M-Wave MAV (No Volitional Trials )",'fontweight','bold','fontsize',14)

%% 
% # Evaluate the filter accuracy




