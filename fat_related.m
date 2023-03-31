%% Fatigue Analysis

%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "mar7" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% Frequency Shift

CondLabels={'Volitional','Stimulated','StimVoli'};
FeatLabels={'MAV','MedFreq','MeanFreq','Ssc','Zc'};
MatLabels={'Unfilt','CombvEMG','CombMwave','GsvEMG','GsMwaves'};
FeatNum=length(FeatLabels);
FeatInd=[1 2 3 4 5];
TrialNum=S.(ExpLabel).iTrial-1;
iVoli=1;
iStim=2;
iStimVoli=3;
AvgLength=0.4;
VoliStart=14;
VoliEnd=VoliStart+AvgLength;
StimStart=29;
StimEnd=StimStart+AvgLength;
StimVoliStart=40;
StimVoliEnd=StimVoliStart+AvgLength;
IndVoli=[VoliStart*StimFreq:VoliEnd*StimFreq];
IndStim=[StimStart*StimFreq:StimEnd*StimFreq];
IndStimVoli=[StimVoliStart*StimFreq:StimVoliEnd*StimFreq];
Ind=[IndVoli;IndStim;IndStimVoli];

for iFeat=1:FeatNum 
    FeatLabel2=FeatLabels{iFeat};
    
    for iFilt=1:FiltNum
        FiltLabel=FiltLabels{iFilt};
%         FeatLabel=sprintf('%svEMGFeatMat',FiltLabel);
        FeatLabel='UnfiltFeatMat';
        
        for iCond=1:3
            CondLabel=CondLabels{iCond};
            
            for iTrial=1:TrialNum
                TrialLabel=sprintf('Trial_%d',iTrial);

                x=S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(Ind(iCond,:),FeatInd(iFeat));
                AvgTrials(iCond,iTrial)= mean(x);
                y=S.(ExpLabel).(TrialLabel).ForceMat(Ind(iCond,:));
                AvgForce(iCond,iTrial)=mean(y);

                avg_x=S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(Ind(iCond,:),iFeat+FeatNum);
                avg_y=S.(ExpLabel).(TrialLabel).AvgForceFrame(Ind(iCond,:),:);

                Mdl = fitlm(avg_y,avg_x);
                Mdlr_sqr=Mdl.Rsquared.Ordinary;
                Mdlcoef_p_val=Mdl.Coefficients.pValue;
                Mdlcoef=Mdl.Coefficients;
                AnovaTable=anova(Mdl,'summary');
                F=table2array(AnovaTable(2,4));
                p_val=table2array(AnovaTable(2,5));
                FeatForceCorr=corr2(avg_y,avg_x);

                S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel2).(CondLabel).Mdlr_sqr=Mdlr_sqr;
                S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel2).(CondLabel).Mdlcoef=Mdlcoef;
                S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel2).(CondLabel).p_val=p_val;
                S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel2).(CondLabel).iCond=iCond;
                S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel2).(CondLabel).FeatForceCorr=FeatForceCorr;

            end
            
            ForMdl = fitlm(AvgForce(iCond,:),AvgTrials(iCond,:));
            ForMdlr_sqr=ForMdl.Rsquared.Ordinary;
            ForMdlcoef_p_val=ForMdl.Coefficients.pValue;
            ForMdlcoef=ForMdl.Coefficients;
            AnovaTable=anova(ForMdl,'summary');
            F=table2array(AnovaTable(2,4));
            Forp_val=table2array(AnovaTable(2,5));
            
            ForceCorr=corr2(AvgForce(iCond,:),AvgTrials(iCond,:));
            
            S.(ExpLabel).(FiltLabel).(FeatLabel2).ForceCorr(iCond)=ForceCorr;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).ForMdlr_sqr(iCond)=ForMdlr_sqr;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).Forp_val(iCond)=Forp_val;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).AvgForce(iCond,:)=AvgForce(iCond,:);
            
            Mdl = fitlm([1:TrialNum],AvgTrials(iCond,:));
            Mdlr_sqr(iCond)=Mdl.Rsquared.Ordinary;
            Mdlcoef_p_val=Mdl.Coefficients.pValue;
            Mdlcoef=Mdl.Coefficients;
            AnovaTable=anova(Mdl,'summary');
            F=table2array(AnovaTable(2,4));
            p_val(iCond)=table2array(AnovaTable(2,5));
            
            PercentShift=(AvgTrials(:,1)-AvgTrials(:,end))./AvgTrials(:,1)*100;
%            MeanShift=(AvgMeanFreqTrials(:,1)-AvgMeanFreqTrials(:,end))./AvgMeanFreqTrials(:,1)*100;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).AvgTrials=AvgTrials;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).AvgForce=AvgForce;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).PercentShift=PercentShift;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).(CondLabel).Mdl=Mdl;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).Mdlr_sqr(iCond)=Mdlr_sqr(iCond);
            S.(ExpLabel).(FiltLabel).(FeatLabel2).(CondLabel).Mdlcoef_p_val=Mdlcoef_p_val;
            S.(ExpLabel).(FiltLabel).(FeatLabel2).p_val(iCond)=p_val(iCond);
            S.(ExpLabel).(FiltLabel).(FeatLabel2).(CondLabel).Mdlcoef=Mdlcoef;
        end
    end
end

%% "Shifts" by fatigue Analysis
ExpLabel='FatigueTrials';
CondLabels={'Volitional', 'Stimulated','StimVoli'};
FeatLabels={'MAV','MedFreq','MeanFreq','Ssc','Zc'};
IndLabels={'Unfilt','CombvEMG','CombMwave','GsvEMG','GsMwave'};
CondLabels={'Volitional', 'Stimulated','StimVoli'};

FeatNum=length(FeatLabels);
FeatInd=[1 2 3 4 5];
TrialNum=S.(ExpLabel).iTrial-1;

% VoliStart=14;
% VoliEnd=14.4;
% StimStart=25;
% StimEnd=25.4;
% StimVoliStart=40;
% StimVoliEnd=40.4;
% IndVoli=[VoliStart*StimFreq:VoliEnd*StimFreq];
% IndStim=[StimStart*StimFreq:StimEnd*StimFreq];
% IndStimVoli=[StimVoliStart*StimFreq:StimVoliEnd*StimFreq];
% IndCond=[IndVoli;IndStim;IndStimVoli];

iCond=2; %For stim-only
iFeat=1; % for MAV
iInd=3; % for *** data
CondLabel=CondLabels{iCond};
FeatLabel2=FeatLabels{iFeat};

FeatIndVec=[1 2 3 4 5 6 7 8 9 10];
iMAV=FeatIndVec(1);
iMedFreq=FeatIndVec(2);
iMeanFreq=FeatIndVec(3);
iSsc=FeatIndVec(4);
iZc=FeatIndVec(5);
iAvgMAV=FeatIndVec(6);
iAvgMedFreq=FeatIndVec(7);
iAvgMeanFreq=FeatIndVec(8);
iAvgSsc=FeatIndVec(9);
iAvgZc=FeatIndVec(10);

ForceBeg=5;
ForceEnd=40;
ForceTimeInd=round([ForceBeg*fs:ForceEnd*fs]);
MovAvgTime=0.5;
MovAvgLength=MovAvgTime*fs;

% StimBeg=25;
% StimEnd=25.4;
% StimTimeInd=[StimBeg*StimFreq:StimEnd*StimFreq];
% VoliBeg=14;
% VoliEnd=14.4;
% VoliTimeInd=[VoliBeg*StimFreq:VoliEnd*StimFreq];

clear AvgMedFreq AvgMAV MeanForceByFrame AvgMeanFreq ttl2 ttl1 ttl
figure(10)
subplot(1,3,3)

for iFilt=1:1
    FiltLabel=FiltLabels{iFilt};
    
    MatLabel=sprintf('%sFrameMat',IndLabels{iInd});
    FeatMatLabel=sprintf('%sFeatMat',IndLabels{iInd});
    FeatLabel=sprintf('%sUnfiltFeatMat',FiltLabel);
%         FeatLabel='UnfiltFeatMat';

    TrialLabel=sprintf('Trial_%d',iTrial);

%         x=S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(IndCond(iCond,:),FeatInd(iFeat));
%         y=S.(ExpLabel).(TrialLabel).ForceMat(IndCond(iCond,:));
%         AvgForce(iCond,iTrial)=mean(y);

%         Force=S.(ExpLabel).(TrialLabel).data(ForceTimeInd,iForce);
%         MAForce=movmean(Force,[MovAvgLength-1 0 ]);
%         S.(ExpLabel).ShiftTrials.MAForce(:,iTrial)=MAForce;

%         AvgMAV(iTrial)=mean(S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(VoliTimeInd,iAvgMAV));
%         AvgMedFreq(iTrial)=mean(S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(VoliTimeInd,iAvgMedFreq));
%         AvgMeanFreq(iTrial)=mean(S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(VoliTimeInd,iAvgMeanFreq));
%         MeanForceByFrame(iTrial)=mean(S.(ExpLabel).(TrialLabel).MeanForceByFrame(StimTimeInd));

    MAV=S.(ExpLabel).(FiltLabel).MAV.AvgTrials(1,:);
    MedFreq=S.(ExpLabel).(FiltLabel).MedFreq.AvgTrials(2,:);
    MeanFreq=S.(ExpLabel).(FiltLabel).MeanFreq.AvgTrials(2,:);
    Force=S.(ExpLabel).(FiltLabel).MeanFreq.AvgForce(2,:);

    r_sqr_MAV=S.(ExpLabel).(FiltLabel).MAV.Mdlr_sqr(1);
    r_sqr_MedFreq=S.(ExpLabel).(FiltLabel).MedFreq.Mdlr_sqr(2);
    r_sqr_MeanFreq=S.(ExpLabel).(FiltLabel).MeanFreq.Mdlr_sqr(2);

    Mdl = fitlm([1:TrialNum],Force);
    r_sqr_Force=Mdl.Rsquared.Ordinary;

    plot([1:length(MedFreq)],MedFreq/MedFreq(1),'LineWidth',2)
    hold on
    plot([1:length(MedFreq)],MeanFreq/MeanFreq(1),'LineWidth',2)
    plot([1:length(MedFreq)],MAV/MAV(1),'LineWidth',2)
    plot([1:length(MedFreq)],Force/Force(1),'LineWidth',2)
    ttl1=sprintf('Voli: %.1f - %.1f secs, Stim: %.1f - %.1f secs',VoliStart,VoliEnd,StimStart,StimEnd);
    ttl2=sprintf('Rsqr: %.2f, %.2f, %.2f, %.2f',r_sqr_MedFreq,r_sqr_MeanFreq,r_sqr_MAV,r_sqr_Force);
    ttl={ttl1; ttl2};
    title( ttl);
    grid

end

legend({'MedFreq (Stim)', 'MeanFreq (Stim)','MAV (Voli)' ,' Force (Stim)' })
xlabel('Fatiguing Trials')
ylabel('Normalized Measures')


%% Fatigue Plotting
PlotTrial=[14 14];
PlotFrame=[175 1400];
FrameInd=[PlotFrame(1):PlotFrame(2)];
FeatPlotLabels={'MAV','Median Freq.','Mean Freq.','SSC','ZC'};
ColorVec={'k','r'};
UnfiltLabel='UnfiltFeatMat';
iExp=4;
ExpLabel=ExpLabels{iExp}; 
TrialNum=S.(ExpLabel).iTrial-1;

for iTrial=PlotTrial(1):PlotTrial(2)
    TrialLabel=sprintf('Trial_%d',iTrial);
    
    for iFeat=1:FeatNum
        FeatPlotLabel=FeatPlotLabels{iFeat};
    
        for iFilt=1:FiltNum
            FiltLabel=FiltLabels{iFilt};
            vEMGLabel=sprintf('%svEMGFeatMat',FiltLabel);
            mWaveLabel=sprintf('%sMwaveFeatMat',FiltLabel);
            
            iLow=FrameInd(1)*FrameLength;
            iHigh=FrameInd(2)*FrameLength;
            
            vEMGFeat=S.(ExpLabel).(TrialLabel).(FiltLabel).(vEMGLabel)(FrameInd,iFeat);
            mWaveFeat=S.(ExpLabel).(TrialLabel).(FiltLabel).(mWaveLabel)(FrameInd,iFeat);
            
            Time=S.(ExpLabel).(TrialLabel).data(iLow:iHigh,iTime);
            figure(1)
            subplot(FiltNum,FeatNum,(iFilt-1)*(FeatNum)+iFeat)
            plot(FrameInd,vEMGFeat,ColorVec{iFilt},'LineWidth',2)
            hold on
            ttl1=sprintf('vEMG %s at Time(s): %.2f-%.2f',FeatPlotLabel,PlotFrame(1)/StimFreq,PlotFrame(2)/StimFreq);
            ttl2=sprintf(' %s, TrialNum: %d',ExpLabel,iTrial);
            title({ttl1; ttl2});
            xlabel('Frames')
            ylabel(FeatPlotLabel)
            legend(FiltLabels{iFilt})
            
            figure(2)
            subplot(FiltNum,FeatNum,(iFilt-1)*(FeatNum)+iFeat)
            plot(FrameInd,mWaveFeat,ColorVec{iFilt},'LineWidth',2)
            hold on
            ttl1=sprintf('m-Wave %s at  Time(s): %.2f-%.2f',FeatPlotLabel,PlotFrame(1)/StimFreq,PlotFrame(2)/StimFreq);
            ttl2=sprintf(' %s, TrialNum: %d',ExpLabel,iTrial);
            title({ttl1;ttl2});
            xlabel('Frames')
            ylabel(FeatPlotLabel)
            legend(FiltLabels{iFilt})
            
        end
        
        Unfilt=S.(ExpLabel).(TrialLabel).(FiltLabel).(UnfiltLabel)(FrameInd,iFeat);
        figure(3)
        subplot(1,FeatNum,iFeat)
        plot(FrameInd,Unfilt,ColorVec{1},'LineWidth',2)
        hold on
        ttl1=sprintf('Unfiltered %s at  Time(s): %.2f-%.2f',FeatPlotLabel,PlotFrame(1)/StimFreq,PlotFrame(2)/StimFreq);
        ttl2=sprintf(' %s, TrialNum: %d',ExpLabel,iTrial);
        title({ttl1;ttl2});
        xlabel('Frames')
        ylabel(FeatPlotLabel)
         
    end                 
end


%%
%Features by trials

EMGLabel={'Mwave','vEMG'};
% FiltPlots={'Unfiltered ','Mwave filtered','vEMG'}
iLbl=1;
for iFeat=1:FeatNum 
    FeatLabel2=FeatLabels{iFeat};
    
    for iFilt=1:FiltNum
        FiltLabel=FiltLabels{iFilt};
        FeatLabel=sprintf('%s%sFeatMat',FiltLabel,EMGLabel{iLbl});
        
        AvgTrials=S.(ExpLabel).(FiltLabel).(FeatLabel2).AvgTrials;
        PercentShift=S.(ExpLabel).(FiltLabel).(FeatLabel2).PercentShift;
        Mdlr_sqr=S.(ExpLabel).(FiltLabel).(FeatLabel2).Mdlr_sqr;
        p_val=S.(ExpLabel).(FiltLabel).(FeatLabel2).p_val;
        
        figure(1)
        subplot(FiltNum,FeatNum,(iFilt-1)*(FeatNum)+iFeat)
        plot(AvgTrials','LineWidth',2)
        xlbl=sprintf('%s filtered %s',FiltLabel,FeatLabel2);
        xlabel('Trials')
        ylabel(xlbl)
%         legend(CondLabels)
        ttl1=sprintf('Trial 1,14 %%Shift:%.2f, %.2f, %.2f',...
            PercentShift(1),PercentShift(2),PercentShift(3));
        ttl2=sprintf('Rsqr:%.2f, %.2f, %.2f, p:%.2f, %.2f, %.2f',...
            Mdlr_sqr(1),Mdlr_sqr(2),Mdlr_sqr(3), p_val(1), p_val(2), p_val(3));
        ttl={ttl1;ttl2};
        title(ttl);
        hold
        
    end
end


subplot(FiltNum,FeatNum,1)
legend(CondLabels)
