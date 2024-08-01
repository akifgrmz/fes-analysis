

clc
clear all
TestFolders=["feb28_24" "feb29_24" "mar18_24" "mar20_24"];
%%
TestFolders=["mar18_24"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% 

clear MAV_Mean
RCMeanTime=[8 10]; % Calculating the means at time [8 10]
FiltLabel="Unfilt";
RampTimeRange=[2 20];

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpRuns=S.(TestStruct).ExpRuns;
    
    
    if ExpRuns(S.(AnaStruct).AnaPar.ExpTable.('RC'))
        ExpLabel=ExpLabels(S.(AnaStruct).AnaPar.ExpTable.('RC'));

        stim_freq=S.(TestStruct).ExpPar.stim_freq;
        MeanFrame=[RCMeanTime(1)*stim_freq RCMeanTime(2)*stim_freq];
        MeanRangeInd=[MeanFrame(1): MeanFrame(2)];

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
        S.(AnaStruct).(ExpLabel).MAV_Vals=table(PWPoints, MAV_Vals(1,:)',...
            MAV_Vals(2,:)',MAV_Vals(3,:)',mean(MAV_Vals)',std(MAV_Vals)',...
            'VariableNames',["PW" "Rep_1" "Rep_2" "Rep_3" "Mean" "Std"]);

        S.(AnaStruct).(ExpLabel).MAV_Mean=array2table(MAV_Mean(:,iTest),...
            'VariableName',TestFolders(iTest));
    end
    

    if ExpRuns(S.(AnaStruct).AnaPar.ExpTable.('Ramp'))
        ExpLabel=S.(AnaStruct).AnaPar.ExpTable.('Ramp');   
        
        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
        MVCLabel=S.(AnaStruct).AnaPar.ExpTable.('MVC');
        OccLabel=S.(AnaStruct).AnaPar.ExpTable.('Occ');

        MVC=S.(TestStruct).(MVCLabel).MVC;
        PercentMVCVals=MVC*S.(TestStruct).(OccLabel).VoliMVCVec/100;
        
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            TrialFrameLength=length(S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames);
            stim_freq=S.(TestStruct).ExpPar.stim_freq;

            FrameRange=stim_freq*RampTimeRange;
            FrameRange(FrameRange>TrialFrameLength)=TrialFrameLength;
            FrameRangeInd=FrameRange(1):FrameRange(2);
            
            MAV=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Filt_MAV_vEMG')(FrameRangeInd);
            Force=S.(AnaStruct).(ExpLabel).(TrialLabel).FiltForceFrames(FrameRangeInd);
            
            for iVal=1:length(PercentMVCVals)
                
                [ValMin(iVal), IndMin(iVal)]=min(abs(PercentMVCVals(iVal)-Force));
                RangeLow=IndMin(iVal)-100;
                RangeHign=IndMin(iVal)+100;
                RangeLow(RangeLow<1)=[];
                RangeHigh(RangeHigh>length(MAV))=[];
                RampMAVVals(:,iVal)=MAV(RangeLow:RangeHigh);
            end
            
            AvgRampMAV=mean(RampMAVVals);
            SDRampMAV=std(RampMAVVals);
            S.(AnaStruct).(ExpLabel).(TrialLabel).RampMAVVals=RampMAVVals;
            S.(AnaStruct).(ExpLabel).(TrialLabel).ValMin=ValMin;
            S.(AnaStruct).(ExpLabel).(TrialLabel).IndMin=IndMin;
            S.(AnaStruct).(ExpLabel).(TrialLabel).AvgRampMAV=AvgRampMAV;
            S.(AnaStruct).(ExpLabel).(TrialLabel).SDRampMAV=SDRampMAV;
            S.(AnaStruct).(ExpLabel).(TrialLabel).PercentMVCVals=PercentMVCVals;
        end
    end 
end



% Update this to more general
% RC_MAV=array2table(MAV_Mean,'VariableNames',TestFolders);
