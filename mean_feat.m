%% Mean MAV of RC Trials 

S = load_test;
%%
clear MAV_Mean
exp_lbl='RC';
MeanTime=[8 10]; % Calculating the means at time [8 10]

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);
    stim_freq=S.(TestStruct).ExpPar.FreqList(1);

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
end

MAV_Mean_Table=array2table(MAV_Mean,'VariableNames',[TestFolders(1),TestFolders(2),TestFolders(3)]);

MAV_Coefs=MAV_Mean./MAV_Mean(end,:);

array2table(MAV_Coefs,'VariableNames',[TestFolders(1),TestFolders(2),TestFolders(3)])
figure
plot(linspace(10,30,7),MAV_Coefs)
legend([TestFolders(1),TestFolders(2),TestFolders(3)])
xlabel('MVC')
ylabel('Norm MAV')

% MAV_Mean=array2table(MAV_Mean,'VariableNames',[sprintf('PW_%d',),TestFolders(2),TestFolders(3)]);
% Normalized MAV values indicated that the normalizing the EMG signals
% among participants are hard to model. 
% But the sample is 3; the third one might an outlier.

%% Mean MAV of Occ Trials 
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
                MAV_Mean(TrialsInd(iTrial,c),2:4)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);
            end
            MAV_Mean_Reps(c,:) = [ c mean(MAV_Mean(TrialsInd(:,c),1)) RepTableMat(TrialsInd(iTrial,c),[4 5 ])];
            c=c+1;
        end
    end
    
    MAV_Mean=array2table(MAV_Mean,'VariableNames',["MAV_Mean","Theo_Force","vMVC","sMVC"]);
    MAV_Mean_Reps=array2table(MAV_Mean_Reps,'VariableNames',["Unique_Trial" "MAV_Mean","vMVC" "sMVC"]);
    
    S.(AnaStruct).(ExpLabel).MAV_Mean=MAV_Mean;
    S.(AnaStruct).(ExpLabel).MAV_Mean_Reps=MAV_Mean_Reps;
    S.(AnaStruct).(ExpLabel).RepTrialsInd=[array2table([TrialsInd'],...
        'VariableNames',["Rep1" "Rep2" "Rep3"]) MAV_Mean_Reps(:,'sMVC') MAV_Mean_Reps(:,'vMVC')];
end
%%
colorcode=[0 0.4470 0.7410; 
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330];

sMVC=0;
for iTest=1:length(TestFolders)
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    sMVCzero_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('sMVC');
    MAVMean_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('MAV_Mean');
    vMVCVal_Reps=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('vMVC');
    
    sMVCzero=S.(AnaStruct).(ExpLabel).MAV_Mean.('sMVC');
    MAVMean=S.(AnaStruct).(ExpLabel).MAV_Mean.('MAV_Mean');
    vMVCVal=S.(AnaStruct).(ExpLabel).MAV_Mean.('vMVC');
    
    Ind=(sMVCzero==sMVC);
    Ind_Reps=(sMVCzero_Reps==sMVC);

    figure(10)
    plot(vMVCVal(Ind),MAVMean(Ind),'*','LineWidth',1,'Color',colorcode(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),'LineWidth',2,'Color',colorcode(iTest,:),'DisplayName',TestFolders(iTest))


end

xlabel('% Voli. MVC')
ylabel('Avg MAV')
title(sprintf('Avg MAV (constant region) with Stim=%d%% MVC',sMVC))
grid
xlim([0 50])