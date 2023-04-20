%% Fatigue Related Analysis
%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);
%% Initial Plotting 
close all
iTest=1;  % pick a test to plot 
FolderName=TestFolders(iTest);  %% Folders to be loaded 
Trials=[1];  % Trial
TimeRange=[5 15];  % in seconds
Exp= ["Fat"];
AnaStruct=sprintf("%s_ana",TestFolders{iTest});
TestStruct=sprintf("%s_test",TestFolders{iTest});
ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(Exp);

NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
if sum(Trials>NumofTrials)>0
    error(sprintf('Trial %d does not exist \n',Trials(Trials>NumofTrials)))
end

for iTrial=Trials
    TrialLabel=sprintf('Trial_%d',iTrial);
    
    Time=S.(AnaStruct).(ExpLabel).(TrialLabel).data.("Time");
    TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
    
    EMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).("BPFilt_EMG");
    Trigger=S.(TestStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).("Trigger");
    Force=S.(TestStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).("Force");
    PW=S.(TestStruct).(ExpLabel).(TrialLabel).data(TimeInd,:).("PW");
    
    figure(iTrial)
    subplot(2,1,1)
    plot(Time(TimeInd),Trigger/1000,'r','LineWidth',1)
    hold
    plot(Time(TimeInd),EMG,'b','LineWidth',2)
    legend({'Trigger(a.u.)','Raw EMG (V)'})
    title('Trigger and EMG Signal')
    xlabel('Time (s)')
    
    subplot(2,1,2)
    plot(Time(TimeInd),PW*100,'r','LineWidth',1)
    hold
    plot(Time(TimeInd),Force,'b','LineWidth',2)
    legend({'PW (a.u.)','Force(N)'})
    title('Pulse Width and Force Signal')
%     set(gca,'XTick',[5 8 11 12 15 18 21])
    xlabel('Time (s)')
end

%% Linear Models
% # Calculating the relevant fitting stats and plotting them
% # Plotting the Mean, Median and MAV 
% # 
Exp= ["Fat"];

VoliStart=5;
VoliEnd=VoliStart+10;
StimStart=15;
StimEnd=StimStart+10;
StimVoliStart=25;
StimVoliEnd=StimVoliStart+10;
Segments=[VoliStart VoliEnd; StimStart StimEnd; StimVoliStart StimVoliEnd];  %% segment beggining and ends 
SegLabels=["Voli" "Stim" "Stim+Voli"];
TestStruct=sprintf("%s_test",TestFolders{1});
AnaStruct=sprintf("%s_ana",TestFolders{1});

stim_freq=S.(TestStruct).ExpPar.stim_freq;
IndVoli=[VoliStart*stim_freq:VoliEnd*stim_freq];
IndStim=[StimStart*stim_freq:StimEnd*stim_freq];
IndStimVoli=[StimVoliStart*stim_freq:StimVoliEnd*stim_freq];
SegInd=[IndVoli;IndStim;IndStimVoli];
[NumofSegs,~]=size(SegInd);

ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(Exp);
FeatLabels=string(S.(AnaStruct).AnaPar.FeatLabels);
FiltLabels=string(S.(AnaStruct).AnaPar.FiltLabels);


Fat_stats=table([],[],[],[],[],[],[],[],[],[],[],[],[],...
    'VariableNames',["Test" "Exp" " Filt" "Trial"...
    "Feat" "Segment" "EMG" "R_sqr" "Coef_p" "Coefs" "Anova_p" "Pear_r" "LPFilt"]);

Fat=table([],[],[],[],[],[],[],[],[],[],[],[],...
    'VariableNames',["Feat_Val" "LPFilt_Val" " Force" "LPFilt_Force"...
    "Frame" "Feat" "Segment" "EMG" "Trial" "Filt" "Test" "Exp"]);
                    
count=1;
for iTest=1:length(TestFolders)
    AnaStruct=sprintf("%s_ana",TestFolders(iTest));
    TestStruct=sprintf("%s_test",TestFolders(iTest));
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
   
    for iFilt=1:length(FiltLabels)
        FiltLabel=FiltLabels(iFilt);            
        samp_filt=1;

        for iTrial=1:NumofTrials
            TrialLabel=sprintf("Trial_%d",iTrial);
            
            for iSeg=1:NumofSegs
                SegLabel=SegLabels(iSeg);
                y=mean(S.(AnaStruct).(ExpLabel).(TrialLabel).ForceFrames(:,SegInd(iSeg,:)));
                y_filt=S.(AnaStruct).(ExpLabel).(TrialLabel).FiltForceFrames(:,SegInd(iSeg,:));
                                
                for iFeat =1:length(FeatLabels)
                    FeatVar=sprintf("%s_vEMG",FeatLabels(iFeat)); 
                
                    x=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats(SegInd(iSeg,:),:).(FeatVar);
                    
                    Mdl = fitlm(y,x); 
                    Mdlr_sqr=Mdl.Rsquared.Ordinary;
                    Mdlcoef_p_val=Mdl.Coefficients.pValue;
                    Mdlcoef=Mdl.Coefficients;
                    AnovaTable=anova(Mdl,'summary');
                    F=table2array(AnovaTable(2,4));
                    p_val=table2array(AnovaTable(2,5));
                    FeatForceCorr=corr2(y',x);

                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlr_sqr(iSeg)=Mdlr_sqr;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef_p_val(:,iSeg)=Mdlcoef_p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef=Mdlcoef;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).p_val(iSeg)=p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).ForceCorr(iSeg)=FeatForceCorr;
                    
                    frames(samp_filt)="Unfilt";
                    emgs(samp_filt)="vEMG";
                    coefs(:,samp_filt)=Mdlcoef.('Estimate');
                    r_sqr(samp_filt)=Mdlr_sqr;
                    coef_p_val(:,samp_filt)=Mdlcoef_p_val;
                    anova_p(samp_filt)=p_val;
                    force_corr(samp_filt)=FeatForceCorr;
                    trials(samp_filt)=TrialLabel;
                    feats(samp_filt)=FeatLabels(iFeat);
                    filts(samp_filt)=FiltLabel;
                    tests(samp_filt)=TestFolders(iTest);
                    exps(samp_filt)=ExpLabel;
                    segs(samp_filt)=SegLabel;
                    samp_filt=samp_filt+1;
                    
                    FeatVar=sprintf("%s_MWave",FeatLabels(iFeat)); 
                    FiltFeatVar=sprintf("Filt_%s_MWave",FeatLabels(iFeat)); 
                    x=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats(SegInd(iSeg,:),:).(FeatVar);
                    
                    Mdl = fitlm(y,x); 
                    Mdlr_sqr=Mdl.Rsquared.Ordinary;
                    Mdlcoef_p_val=Mdl.Coefficients.pValue;
                    Mdlcoef=Mdl.Coefficients;
                    AnovaTable=anova(Mdl,'summary');
                    F=table2array(AnovaTable(2,4));
                    p_val=table2array(AnovaTable(2,5));
                    FeatForceCorr=corr2(y',x);

                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlr_sqr(iSeg)=Mdlr_sqr;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef_p_val(:,iSeg)=Mdlcoef_p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef=Mdlcoef;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).p_val(iSeg)=p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).ForceCorr(iSeg)=FeatForceCorr;
                    
                    frames(samp_filt)="Unfilt";
                    emgs(samp_filt)="MWave";
                    coefs(:,samp_filt)=Mdlcoef.('Estimate');
                    r_sqr(samp_filt)=Mdlr_sqr;
                    coef_p_val(:,samp_filt)=Mdlcoef_p_val;
                    anova_p(samp_filt)=p_val;
                    force_corr(samp_filt)=FeatForceCorr;
                    trials(samp_filt)=TrialLabel;
                    feats(samp_filt)=FeatLabels(iFeat);
                    filts(samp_filt)=FiltLabel;
                    tests(samp_filt)=TestFolders(iTest);
                    exps(samp_filt)=ExpLabel;
                    segs(samp_filt)=SegLabel;
                    samp_filt=samp_filt+1;
                    
                    FiltFeatVar=sprintf("Filt_%s_vEMG",FeatLabels(iFeat));
                    FeatVar=sprintf("%s_vEMG",FeatLabels(iFeat)); 
                    x_filt=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(SegInd(iSeg,:),:).(FiltFeatVar);

                    Mdl = fitlm(y_filt,x_filt); 
                    Mdlr_sqr=Mdl.Rsquared.Ordinary;
                    Mdlcoef_p_val=Mdl.Coefficients.pValue;
                    Mdlcoef=Mdl.Coefficients;
                    AnovaTable=anova(Mdl,'summary');
                    F=table2array(AnovaTable(2,4));
                    p_val=table2array(AnovaTable(2,5));
                    FeatForceCorr=corr2(y',x);

                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlr_sqr(iSeg)=Mdlr_sqr;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef_p_val(:,iSeg)=Mdlcoef_p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef=Mdlcoef;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).p_val(iSeg)=p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).ForceCorr(iSeg)=FeatForceCorr;
                    
                    frames(samp_filt)="LPFilt";
                    emgs(samp_filt)="vEMG";
                    coefs(:,samp_filt)=Mdlcoef.('Estimate');
                    r_sqr(samp_filt)=Mdlr_sqr;
                    coef_p_val(:,samp_filt)=Mdlcoef_p_val;
                    anova_p(samp_filt)=p_val;
                    force_corr(samp_filt)=FeatForceCorr;
                    trials(samp_filt)=TrialLabel;
                    feats(samp_filt)=FeatLabels(iFeat);
                    filts(samp_filt)=FiltLabel;
                    tests(samp_filt)=TestFolders(iTest);
                    exps(samp_filt)=ExpLabel;
                    segs(samp_filt)=SegLabel;
                    samp_filt=samp_filt+1;
                    
                    FiltFeatVar=sprintf("Filt_%s_MWave",FeatLabels(iFeat)); 
                    FeatVar=sprintf("%s_MWave",FeatLabels(iFeat)); 
                    x_filt=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(SegInd(iSeg,:),:).(FiltFeatVar);
                    
                    Mdl = fitlm(y_filt,x_filt); 
                    Mdlr_sqr=Mdl.Rsquared.Ordinary;
                    Mdlcoef_p_val=Mdl.Coefficients.pValue;
                    Mdlcoef=Mdl.Coefficients;
                    AnovaTable=anova(Mdl,'summary');
                    F=table2array(AnovaTable(2,4));
                    p_val=table2array(AnovaTable(2,5));
                    FeatForceCorr=corr2(y',x);

                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlr_sqr(iSeg)=Mdlr_sqr;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef_p_val(:,iSeg)=Mdlcoef_p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Mdlcoef=Mdlcoef;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).p_val(iSeg)=p_val;
                    S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).ForceCorr(iSeg)=FeatForceCorr;
                    
                    frames(samp_filt)="LPFilt";
                    emgs(samp_filt)="MWave";
                    coefs(:,samp_filt)=Mdlcoef.('Estimate');
                    r_sqr(samp_filt)=Mdlr_sqr;
                    coef_p_val(:,samp_filt)=Mdlcoef_p_val;
                    anova_p(samp_filt)=p_val;
                    force_corr(samp_filt)=FeatForceCorr;
                    trials(samp_filt)=TrialLabel;
                    feats(samp_filt)=FeatLabels(iFeat);
                    filts(samp_filt)=FiltLabel;
                    tests(samp_filt)=TestFolders(iTest);
                    exps(samp_filt)=ExpLabel;
                    segs(samp_filt)=SegLabel;
                    samp_filt=samp_filt+1;
                        
                    
                    FeatVar=sprintf("%s_vEMG",FeatLabels(iFeat)); 
                    FiltFeatVar=sprintf("Filt_%s_vEMG",FeatLabels(iFeat)); 
                    x=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats(SegInd(iSeg,:),:).(FeatVar);
                    x_filt=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(SegInd(iSeg,:),:).(FiltFeatVar);
                    
                    lg=length(x);
                    g_frame=strings(1,lg);
                    g_frame(:)=1:lg;
                    g_trial=strings(1,lg);
                    g_trial(:)=TrialLabel;
                    g_feat=strings(1,lg);
                    g_feat(:)=FeatLabels(iFeat); 
                    g_filt=strings(1,lg);
                    g_filt(:)=FiltLabel;
                    g_test=strings(1,lg);
                    g_test(:)=TestFolders(iTest);
                    g_exp=strings(1,lg);
                    g_exp(:)=ExpLabel;
                    g_seg=strings(1,lg);
                    g_seg(:)=SegLabel;
                    g_emg=strings(1,lg);
                    g_emg(:)="vEMG";
                    
                    temp=table(x, x_filt, y', y_filt' ,g_frame' ,g_feat' ,g_seg' ,g_emg'...
                        ,g_trial' ,g_filt' ,g_test' ,g_exp','VariableNames',...
                        ["Feat_Val" "LPFilt_Val" " Force" "LPFilt_Force"...
                        "Frame" "Feat" "Segment" "EMG" "Trial" "Filt" "Test" "Exp"]);
                        
                    Fat= [Fat; temp];
                    
                    FeatVar=sprintf("%s_MWave",FeatLabels(iFeat)); 
                    FiltFeatVar=sprintf("Filt_%s_MWave",FeatLabels(iFeat)); 
                    x=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats(SegInd(iSeg,:),:).(FeatVar);
                    x_filt=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(SegInd(iSeg,:),:).(FiltFeatVar);
                    
                    lg=length(x);
                    
                    g_frame=strings(1,lg);
                    g_frame(:)=1:lg;
                    g_trial=strings(1,lg);
                    g_trial(:)=TrialLabel;
                    g_feat=strings(1,lg);
                    g_feat(:)=FeatLabels(iFeat); 
                    g_filt=strings(1,lg);
                    g_filt(:)=FiltLabel;
                    g_test=strings(1,lg);
                    g_test(:)=TestFolders(iTest);
                    g_exp=strings(1,lg);
                    g_exp(:)=ExpLabel;
                    g_seg=strings(1,lg);
                    g_seg(:)=SegLabel;
                    g_emg=strings(1,lg);
                    g_emg(:)="MWave";
                    
                    temp=table(x, x_filt, y', y_filt' ,g_frame' ,g_feat' ,g_seg' ,g_emg'...
                        ,g_trial' ,g_filt' ,g_test' ,g_exp','VariableNames',...
                        ["Feat_Val" "LPFilt_Val" " Force" "LPFilt_Force"...
                        "Frame" "Feat" "Segment" "EMG" "Trial" "Filt" "Test" "Exp"]);
                        
                    Fat= [Fat; temp];

                end
            end
        end

        FatTable=table(tests',exps',filts',trials',feats',segs',emgs',r_sqr',...
            coef_p_val',coefs',anova_p',force_corr',frames',...
            'VariableNames',["Test" "Exp" " Filt" "Trial"...
            "Feat" "Segment" "EMG" "R_sqr" "Coef_p" "Coefs" "Anova_p" "Pear_r" "LPFilt"]);
        S.(AnaStruct).(ExpLabel).(FiltLabel).FatTable=FatTable;
        Fat_stats= [ Fat_stats; FatTable];
        clear tests exps filts trials feats r_sqr coef_p_val anova_p force_corr segs coefs emgs frames FatTable
    
    end
end

writetable( Fat_stats, 'fatigue_stats.csv')
writetable( Fat, 'fatigue_frames.csv')

%% Fatigue Features Shifts
% # Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% Plotting Features  
close all
iTest=1;  % pick a test to plot 
FolderName='jan7';  %% Folders to be loaded 
Trials=[1 5];  % Trial
TimeRange=[1 35];  % in seconds
Exp= ["Fat"];
FiltLabel="Unfilt";
AnaLabel=sprintf("%s_ana",TestFolders{iTest});
TestLabel=sprintf("%s_test",TestFolders{iTest});
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(Exp);
stim_freq=S.(TestLabel).ExpPar.stim_freq;
FrameInd=TimeRange(1)*stim_freq:TimeRange(2)*stim_freq;
cm = lines(length(TestFolders));

NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
if sum(Trials>NumofTrials)>0
    error(sprintf('Trial %d does not exist \n',Trials(Trials>NumofTrials)))
end

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});

    for iTrial=Trials
        TrialLabel=sprintf('Trial_%d',iTrial);

        Frames=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(:,FrameInd));
        Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames(:,FrameInd));
        
        FatMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).("Filt_MAV_vEMG");
        FatMedFreq=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).("Filt_MedFreq_vEMG");
        FatMeanFreq=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).("Filt_MeanFreq_vEMG");

        MAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameInd,:).("MAV_vEMG");
        MedFreq=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameInd,:).("MedFreq_vEMG");
        MeanFreq=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameInd,:).("MeanFreq_vEMG");

        figure(iTrial)
        
        subplot(4,1,1)
        plot(Frames,Force,'LineWidth',2.5,'DisplayName',sprintf("%s Force(N)",TestFolders{iTest}),'Color',cm(iTest,:))
        hold on
        title(sprintf('Force of Trial %d of Different Subjects',iTrial))
        xlabel('Time (s)')
        legend
        grid on
        
        subplot(4,1,2)
        plot(Frames,MAV,'--','LineWidth',0.02,'DisplayName',TestFolders{iTest},'Color',cm(iTest,:))
        hold on 
        plot(Frames,FatMAV,'LineWidth',2.5,'DisplayName',sprintf("%s Filtered",TestFolders{iTest}),'Color',cm(iTest,:))
        title(sprintf('MAV of Trial %d of Different Subjects',iTrial))
        xlabel('Time (s)')
        ylim([0 max(MAV)/3])
        legend
        grid on

        subplot(4,1,3)
        plot(Frames,MedFreq,'--','LineWidth',.02,'DisplayName',TestFolders{iTest},'Color',cm(iTest,:))
        hold on
        plot(Frames,FatMedFreq,'LineWidth',2.5,'DisplayName',sprintf("%s Filtered",TestFolders{iTest}),'Color',cm(iTest,:))
        title(sprintf('Median Freq. of Trial %d of Different Subjects',iTrial))
        xlabel('Time (s)')
        legend
        grid on
        
        subplot(4,1,4)
        plot(Frames,MeanFreq,'--','LineWidth',.02,'DisplayName',TestFolders{iTest},'Color',cm(iTest,:))
        hold on
        plot(Frames,FatMeanFreq,'LineWidth',2.5,'DisplayName',sprintf("%s Filtered",TestFolders{iTest}),'Color',cm(iTest,:))
        title(sprintf('Mean Freq. of Trial %d of Different Subjects',iTrial))
        xlabel('Time (s)')
        ylabel('Hz')
        legend
        grid on 

    end
end

%% # Fatigue by Trials Plotting 
clc
Exp="Fat";
SegLabels=["Voli" "Stim" "Stim+Voli"];

FiltLabel="Unfilt";
EMGLabel="vEMG";
FeatFilt="FatFeats";
FeatFiltLabel="Filt_";

VoliStart=10;
VoliEnd=VoliStart+5;
StimStart=20;
StimEnd=StimStart+5;
StimVoliStart=30;
StimVoliEnd=StimVoliStart+5;

AnaLabel=sprintf("%s_ana",TestFolders(1));
TestLabel=sprintf("%s_test",TestFolders{1});
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(Exp);

stim_freq=S.(TestLabel).ExpPar.stim_freq;
IndVoli=[VoliStart*stim_freq:VoliEnd*stim_freq];
IndStim=[StimStart*stim_freq:StimEnd*stim_freq];
IndStimVoli=[StimVoliStart*stim_freq:StimVoliEnd*stim_freq];
SegInd=[IndVoli;IndStim;IndStimVoli];
[NumofSegs,~]=size(SegInd);

for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders(iTest));
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    
    FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
    FeatLabels=string(S.(AnaLabel).AnaPar.FeatLabels);
    NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
    clear FiltFeat
    
    FatTable=S.(AnaLabel).(ExpLabel).(FiltLabel).FatTable;

    [SegNum,~]=size(SegInd);
    for iSeg=1:SegNum
        SegLabel=SegLabels(iSeg);

        for iFeat=1:3 % MAV, Med, Median Freq
            FeatLabel=sprintf("%s%s_%s",FeatFiltLabel, FeatLabels(iFeat), EMGLabel);

            for iTrial=1:NumofTrials
                TrialLabel=sprintf("Trial_%d",iTrial);
%                     FiltForce=S.(AnaLabel).(ExpLabel).(TrialLabel).FiltForceFrames(SegInd(:,iSeg));
                FiltFeat(iTrial)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(FeatFilt)(SegInd(iSeg,:),:).(FeatLabel));
            end

            figure(100)
            subplot(SegNum,length(TestFolders),iSeg*SegNum-SegNum+iTest)
            plot(1:NumofTrials,FiltFeat/FiltFeat(1),'DisplayName',FeatLabel)
            hold on 
            title(sprintf("Test: %s, Segment: %s, Filt: %s ",TestFolders(iTest),SegLabel,FiltLabel))
            xlabel('Trials')
            ylabel('Norm. EMG Feature')
        end
    end
end

legend('Location', 'NorthWest')

% # Figures with mWaves 

% # Figures with vEMG 

%% "Shifts" by fatigue Analysis

for iFilt=1:1
    FiltLabel=FiltLabels{iFilt};

    TrialLabel=sprintf('Trial_%d',iTrial);

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
