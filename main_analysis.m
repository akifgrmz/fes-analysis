%% mainscript v0.2
% Strcut S contains experimental data without any processing 
% Struct P contains processed data
% this file is for running the analysis for one file only but will be a
% function for running all the files in future
%% Run tidy_data if you have not yet done so 
clear all
FolderName="jan7";
S=tidy_data(FolderName);

%% Data Inject 

clc
clear all
TestFolders=["jan7"];
TestFiles=["jan7_test"];
StructstoLoad=["ExpPar"]; 
TestStruct=TestFiles(1);
AnaStruct=sprintf("%s_ana",TestFolders);
S = load_test(TestFolders,TestFiles);

fs=S.(TestStruct).ExpPar.fs;

%% Some Experimental Parameters Required for the Analysis
% Experimental Parameters



%% Defining the parameters 
% Experiment Parameters

BlankingTime= 0.004; % secs
StimTime=0.007;
TrigDelay=0.006;
TrigTime=0.003;

% Experiments indices 
ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
DataInd= S.(TestStruct).ExpPar.DataInd;
TableInd=DataInd.Properties.VariableNames;
iForce=table2array(DataInd(:,"Force"));
iTrigger=table2array(DataInd(:,"Trigger"));
iEMG=table2array(DataInd(:,"EMG"));
iTime=table2array(DataInd(:,"Time"));
iPW=table2array(DataInd(:,"PW"));

% Analysis Parameters
S.(AnaStruct).AnaPar.AnaLabels=["BPFilt_EMG"];
S.(AnaStruct).AnaPar.AnaInd=table([],'VariableNames',S.(AnaStruct).AnaPar.AnaLabels);

GsOrder=6;

TrigLowThres=0.5;
TrigHighThres=4;
TrigThres=1;

% creating the structs

for iExp=1:length(ExpLabels)
    ExpLabel=ExpLabels{iExp};
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
        S.(AnaStruct).(ExpLabel).(TrialLabel)=struct;
%         S.(AnaStruct).(ExpLabel).(TrialLabel).data=table([],[],[],[],[],'VariableNames',TableInd);

    end
end

%% Design 10-500Hz 20th order butterworth for filtfilt

BPOrder=20;
fcutlow = 10;
fcuthigh = 500;

S.(AnaStruct).AnaPar.BPFilter.BPOrder=BPOrder;
S.(AnaStruct).AnaPar.BPFilter.fcutlow=fcutlow;
S.(AnaStruct).AnaPar.BPFilter.fcuthigh=fcuthigh;

d1 = designfilt('bandpassiir','FilterOrder',BPOrder, ...
         'HalfPowerFrequency1',fcutlow,'HalfPowerFrequency2',fcuthigh, ...
         'SampleRate',fs,'DesignMethod',"butter");

% d2=butter(BPOrder, [fcutlow fcuthigh]/(fs/2), 'bandpass');
% fvtool(d1,'MagnitudeDisplay','magnitude')
%
AnaLabels=S.(AnaStruct).AnaPar.AnaLabels;
AmpGain=990;
for iExp=1:length(ExpLabels)
    ExpLabel=ExpLabels{iExp};
    
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
    
    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
% 
        x=S.(TestStruct).(ExpLabel).(TrialLabel).data.(TableInd{iEMG})/AmpGain;
%         
        BPFilt_EMG = filtfilt(d1,x);
        S.(AnaStruct).(ExpLabel).(TrialLabel).data=table(BPFilt_EMG);
    end
end


%% Trigger 

% P.ExpTrials.RepTableMat=S.ExpTrials.RepTableMat;
BlankTime=0.004;
BlankLength=round(BlankTime*fs);
AnaLabels=S.(TestStruct).ExpPar.ExpLabels
ExpstoAna=AnaLabels([2,4,5]);

for iExp=1:length(ExpstoAna)
    ExpLabel=ExpstoAna{iExp}
 
    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
    
    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        y=S.(AnaStruct).(ExpLabel).(TrialLabel).data.("BPFilt_EMG");
        x=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Trigger");
        
        Time=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Time"); 
        length(y)
            % # blanking
        RisingInd=diff(x)>0;
        FallingInd=diff(x)<0;
        BegofFrames=find(FallingInd>0);
        FallingInd(end+1)=0;
        RisingInd(end+1)=0;

        for iFall=1:length(BegofFrames)-1
            
            y(BegofFrames(iFall):BegofFrames(iFall)+BlankLength)=0;
        end
%         

        S.(AnaStruct).(ExpLabel).(TrialLabel).Indices=table([FallingInd], [RisingInd]);
        S.(AnaStruct).(ExpLabel).(TrialLabel).data.("BlankedEMG")=y;
        S.(AnaStruct).(ExpLabel).(TrialLabel).data.("Time")=Time;


    end
end


%% Plotting after blanking
PlotExp=[ 3 3] ; % 1-'MVCTrials' 2-'RCCurveTrials', 3-'ExpTrials', 4-'FatigueTrials'
PlotTrial=[5 5];
PlotFrame=floor([ 336 337]);
% PlotFrame=floor([ 194 195]);

PlotIndicetype=1; % 1-RawEMG 2-Trigger 3-Force 4-PW 5-Time

for iExp=PlotExp(1):PlotExp(2)
    ExpLabel=ExpLabels{iExp};
    
    TrialNum=S.(ExpLabel).iTrial-1;
    
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        EMGIndiceLow=BlankLength+PlotFrame(1)*FrameLength;
        EMGIndiceHigh=BlankLength+(PlotFrame(2)+1)*FrameLength;
%         FrameNum=floor(length(y)/FrameLength);

        iBPFilterEMG=P.(ExpLabel).(TrialLabel).data(EMGIndiceLow:EMGIndiceHigh,iBPFilter);
        EMG=P.(ExpLabel).(TrialLabel).data(EMGIndiceLow:EMGIndiceHigh,iEMG);
        BlankEMG=P.(ExpLabel).(TrialLabel).data(EMGIndiceLow:EMGIndiceHigh,iBlankEMG);
        Trig=P.(ExpLabel).(TrialLabel).data(EMGIndiceLow:EMGIndiceHigh,iTrig);

        Time=P.(ExpLabel).(TrialLabel).data(EMGIndiceLow:EMGIndiceHigh,iTime);    
        figure
        subplot(2,1,1)
        plot(Time,EMG,'k','LineWidth',2)
        hold on
        plot(Time,BlankEMG,'r','LineWidth',2)
%         plot(Time,Trig/1000,'b','LineWidth',2)

%         plot(Time,Trig,'b','LineWidth',2)

        legend({'Pre-Blanking,BPfiltering', 'Post-Blanking,BPfiltering'})
        ttl=sprintf('BP Filtered %s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{iEMG},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2));
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        ylim([1.5*min(iBPFilterEMG) 1.5*max(iBPFilterEMG)])
        subplot(2,1,2)
        plot(Time,Trig/1000,'b','LineWidth',2)
    end
end


%% Frames Matrix for Force and EMG
% (also remove blanked periods)
% (Add unfilt as a filter and add e m-wave as zero matrix)
% Identify Zero stim trials -> Done

% Identifiying PW=0 frames
for iExp=1:ExpNum
    ExpLabel=ExpLabels{iExp};

    TrialNum=P.(ExpLabel).TrialNum;
    for iTrial=1:TrialNum
        TrialLabel=sprintf('Trial_%d',iTrial);

        PW=P.(ExpLabel).(TrialLabel).data(:,iPW);
        FrameInd=P.(ExpLabel).(TrialLabel).FramesInd;
        FrameNum=length(FrameInd);

        clear PWFrames
        for iFrame=1:FrameNum-1
            %PW
            PWFrames(iFrame)=PW(FrameInd(iFrame));
        end

        P.(ExpLabel).(TrialLabel).PWFrames=PWFrames;
    end
end

%
BlankedFrameLength=79; % extends blanking period to eliminate the trigger variance determined based on trial and error
ForceFrameLength=93;
for iExp=1:ExpNum
    ExpLabel=ExpLabels{iExp};

    TrialNum=P.(ExpLabel).TrialNum;

    for iTrial=1:TrialNum
        TrialLabel=sprintf('Trial_%d',iTrial);

        x=P.(ExpLabel).(TrialLabel).data(:,iBlankEMG);
        z=P.(ExpLabel).(TrialLabel).data(:,iForce);
        PW=P.(ExpLabel).(TrialLabel).data(:,iPW);
        Time=P.(ExpLabel).(TrialLabel).data(:,iTime);

        FrameInd=P.(ExpLabel).(TrialLabel).FramesInd;
        FrameNum=length(FrameInd);
        
        clear y t FrameLengths PWFrames
        for iFrame=1:FrameNum-1
            FrameLengths(iFrame)=-FrameInd(iFrame)+FrameInd(iFrame+1);
            % # EMG
            y(:,iFrame)=x(FrameInd(iFrame+1)-BlankedFrameLength:FrameInd(iFrame+1));
            % # Force
            t(:,iFrame)= z(FrameInd(iFrame+1)-ForceFrameLength:FrameInd(iFrame+1));
            % # PW
%             PWFrames(iFrame)=PW(FrameInd(iFrame));
            % # Time
            TmFrames(:,iFrame)= Time(FrameInd(iFrame+1)-ForceFrameLength:FrameInd(iFrame+1));
        end

        P.(ExpLabel).(TrialLabel).EMGFrames=y;
        P.(ExpLabel).(TrialLabel).ForceFrames=t;
        P.(ExpLabel).(TrialLabel).TimeFrames=TmFrames;
        P.(ExpLabel).(TrialLabel).FrameLengths=FrameLengths;

    end
end

%% Seperating dropped frames

for iExp=1:ExpNum
    ExpLabel=P.Labels.ExpLabels{iExp};
    
    TrialNum=P.(ExpLabel).TrialNum;

    if iExp==3
        StimMVC=P.(ExpLabel).RepTableMat(:,4);
        ZeroStimTrialsInd=find((StimMVC==0));
        DroppedRange1=[5.1 10.8]*StimFreq;
        DroppedRange2=[16.1 18.9]*StimFreq;
        for iTrial=ZeroStimTrialsInd'
            TrialLabel=sprintf('Trial_%d',iTrial);
            
            DroppedFramesEMG=[];
            DroppedFramesForce=[];
            DroppedFramesIndex=[];
            DroppedFramesIndexFixed=DroppedFramesIndex;
            DroppedPWInd=[];
            P.(ExpLabel).(TrialLabel).DroppedFramesEMG=DroppedFramesEMG;
            P.(ExpLabel).(TrialLabel).DroppedFramesForce=DroppedFramesForce;
            P.(ExpLabel).(TrialLabel).DroppedFramesIndex=DroppedFramesIndex;
            P.(ExpLabel).(TrialLabel).DroppedFramesIndexFixed=DroppedFramesIndexFixed;

            P.(ExpLabel).(TrialLabel).DroppedPWInd=DroppedPWInd;        
        end
            
        P.(ExpLabel).ZeroStimTrialsInd=ZeroStimTrialsInd;
    else
        DroppedRange1=[1000 1001];
        DroppedRange2=[1000 1001];
        ZeroStimTrialsInd=[];
    end
    
    for iTrial=setdiff(1:TrialNum,ZeroStimTrialsInd)
        TrialLabel=sprintf('Trial_%d',iTrial);

        Time=P.(ExpLabel).(TrialLabel).data(:,iTime);
        PWFrames=P.(ExpLabel).(TrialLabel).PWFrames;
        x1=P.(ExpLabel).(TrialLabel).EMGFrames;
        x2=P.(ExpLabel).(TrialLabel).ForceFrames;

        FrameNum=length(PWFrames);
        DroppedPWInd1= (DroppedRange1(1)< 1:FrameNum & 1:FrameNum < DroppedRange1(2) & PWFrames == 0 ); 
        DroppedPWInd2= (DroppedRange2(1)< 1:FrameNum & 1:FrameNum < DroppedRange2(2) & PWFrames == 0 ); 
        DroppedPWInd=DroppedPWInd1+DroppedPWInd2;
        DroppedFramesIndex=find(DroppedPWInd==1);
        
        if isempty(DroppedFramesIndex)
            DroppedFramesIndexFixed=DroppedFramesIndex;
        else
            DroppedFramesIndexFixed=[194,230,266,301,336,373,587,623,658];
        end
        
        DroppedFramesEMG=x1(:,DroppedFramesIndexFixed);
        DroppedFramesForce=x2(:,DroppedFramesIndexFixed);

        P.(ExpLabel).(TrialLabel).DroppedFramesEMG=DroppedFramesEMG;
        P.(ExpLabel).(TrialLabel).DroppedFramesForce=DroppedFramesForce;
        P.(ExpLabel).(TrialLabel).DroppedFramesIndex=DroppedFramesIndex;
        P.(ExpLabel).(TrialLabel).DroppedFramesIndexFixed=DroppedFramesIndexFixed;

        P.(ExpLabel).(TrialLabel).DroppedPWInd=DroppedPWInd;

    end
end
%% Test
iExp=4;
ExpLabel=ExpLabels{iExp};
FrameInd=P.(ExpLabel).(TrialLabel).FramesInd;
FrameNum=length(FrameInd);
x1=P.(ExpLabel).(TrialLabel).EMGFrames;
DroppedFramesIndex=[];
y1=x1(:,setdiff(1:FrameNum-1,DroppedFramesIndex));
%% M-wave filtering  
%( Remove blanked periods before filtering, Do unfiltered Frames)
% (update to recreate the vector forms of the signals)

for iExp=1:ExpNum
    ExpLabel=P.Labels.ExpLabels{iExp};
    
    TrialNum=P.(ExpLabel).TrialNum;
    
    for iTrial=1:TrialNum
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        FrameInd=P.(ExpLabel).(TrialLabel).FramesInd;
        FrameNum=length(FrameInd);
        DroppedFramesIndexFixed=P.(ExpLabel).(TrialLabel).DroppedFramesIndexFixed;
        x_frames=P.(ExpLabel).(TrialLabel).EMGFrames(:,setdiff(1:FrameNum-1,DroppedFramesIndexFixed)); % excluding dropped frames 
        [FrameLengthBlanked, FrameNum]=size(x_frames);

        % # Unfilt
        FiltLabel=P.Labels.FiltLabels{3};
        UnfiltvEMGFrames=x_frames;
        UnfiltMWaveFrames=zeros(size(UnfiltvEMGFrames));
        UnfiltvEMG=reshape(UnfiltvEMGFrames,1,FrameLengthBlanked*FrameNum);
        UnfiltMWave=reshape(UnfiltMWaveFrames,1,FrameLengthBlanked*FrameNum);
        
        P.(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=UnfiltMWaveFrames;
        P.(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=UnfiltvEMGFrames;
        P.(ExpLabel).(TrialLabel).(FiltLabel).vEMG=UnfiltvEMG;
        P.(ExpLabel).(TrialLabel).(FiltLabel).MWave=UnfiltMWave;
        
        % # Comb Filter
        [CombMWaveFrames,CombvEMGFrames]=FiltComb(x_frames);
        CombvEMG=reshape(CombvEMGFrames,1,FrameLengthBlanked*FrameNum);
        CombMWave=reshape(CombMWaveFrames,1,FrameLengthBlanked*FrameNum);
        
        FiltLabel=P.Labels.FiltLabels{1};
        P.(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=CombMWaveFrames;
        P.(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=CombvEMGFrames;
        P.(ExpLabel).(TrialLabel).(FiltLabel).vEMG=CombvEMG;
        P.(ExpLabel).(TrialLabel).(FiltLabel).MWave=CombMWave;
        
        % # GS Filter
        % This method is taking too long to calculate
        [~,FrameNum]=size(x_frames);
        clear vEMGMat MWaveMat
        for iFrame = 1:FrameNum-GsOrder
            temp_vect = SUB_GS_filter(x_frames(:,iFrame:GsOrder+iFrame),FrameLength,GsOrder);
            vEMGMat(:,iFrame+GsOrder) = temp_vect(:,1);
            MWaveMat(:,iFrame+GsOrder) = temp_vect(:,2);
        end
        
        FiltLabel=P.Labels.FiltLabels{2};

        P.(ExpLabel).(TrialLabel).(P.Labels.FiltLabels{2}).vEMGFrames=vEMGMat;
        P.(ExpLabel).(TrialLabel).(P.Labels.FiltLabels{2}).MWaveFrames=MWaveMat;
        
        GSvEMG=reshape(vEMGMat,1,FrameLengthBlanked*FrameNum);
        GSMWave=reshape(MWaveMat,1,FrameLengthBlanked*FrameNum);
        P.(ExpLabel).(TrialLabel).(FiltLabel).vEMG=GSvEMG;
        P.(ExpLabel).(TrialLabel).(FiltLabel).MWave=GSMWave;

        % # Blanking filter comes here
                      
    end
end
%% Plotting, Debugging 
PlotExp=[ 3 3] ; % 1-'MVCTrials' 2-'RCCurveTrials', 3-'ExpTrials', 4-'FatigueTrials'
PlotTrial=[5 5];
PlotFrame=floor([ StimFreq*6.8+10 StimFreq*6.8+11]);
PlotFrame=floor([ 336 337]);

PlotIndicetype=1; % 1-RawEMG 2-Trigger 3-Force 4-Time 5-BlankedEMG
FixedFrameLength=80;
for iExp=PlotExp(1):PlotExp(2)
    ExpLabel=P.Labels.ExpLabels{iExp};
    
    TrialNum=P.(ExpLabel).TrialNum;
    
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        
%         FrameLengths=P.(ExpLabel).(TrialLabel).FrameLengths;
%         
%         EMGIndiceLow=BlankLength+PlotFrame(1)*FrameLengths;
%         EMGIndiceHigh=BlankLength+(PlotFrame(2)+1)*FrameLength;
%         FrameNum=floor(length(y)/FrameLength);

        iFilt=1;
        FiltLabel=P.Labels.FiltLabels{iFilt};
        CombvEMGFrames=P.(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        CombMWaveFrames=P.(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        CombvEMG=reshape(CombvEMGFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));
        CombMWave=reshape(CombMWaveFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));

        iFilt=2;
        FiltLabel=P.Labels.FiltLabels{iFilt};
        GSvEMGFrames=P.(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        GSMWaveFrames=P.(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        GSvEMG=reshape(GSvEMGFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));
        GSMWave=reshape(GSMWaveFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));

        UnfiltEMGFrames=P.(ExpLabel).(TrialLabel).EMGFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        UnfiltEMG=reshape(UnfiltEMGFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));

        TimeFrames=P.(ExpLabel).(TrialLabel).TimeFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
        Time=reshape(TimeFrames,1,FixedFrameLength*length(PlotFrame(1):PlotFrame(2)));

        figure
        subplot(2,1,1)
        plot(Time,UnfiltEMG,'k','LineWidth',2)
        hold
        plot(Time,CombvEMG,'r','LineWidth',2)
        legend({'Pre-Comb Filtered', 'Post-Comb Filtered'})
        ttl=sprintf('%s at %s, TrialNum: %d, Frame: %d-%d',P.Labels.DataLabels{8},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        subplot(2,1,2)
        plot(Time,UnfiltEMG,'k','LineWidth',2)
        hold
        plot(Time,CombMWave,'r','LineWidth',2)
        legend({'Pre-Comb Filtered', 'Post-Comb Filtered'})
        ttl=sprintf(' %s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{9},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})

        figure
        subplot(2,1,1)
        plot(Time,UnfiltEMG,'k','LineWidth',2)
        hold
        plot(Time,GSvEMG,'r','LineWidth',2)
        legend({'Pre-GS Filtered', 'Post-GS Filtered'})
        ttl=sprintf('%s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{10},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        subplot(2,1,2)
        plot(Time,UnfiltEMG,'k','LineWidth',2)
        hold
        plot(Time,GSMWave,'r','LineWidth',2)
        legend({'Pre-GS Filtered', 'Post-GS Filtered'})
        ttl=sprintf(' %s at %s, TrialNum: %d, Frame: %d-%d',DataLabels{11},ExpLabel,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        
    end
end

%% EMG Features
% (move force averaging to the top where BP is held)
FeatIndVec=[1 2 3 4 5 6 7 8 9 10];
iMAV=FeatIndVec(1);
iMedFreq=FeatIndVec(2);
iMeanFreq=FeatIndVec(3);
iSsc=FeatIndVec(4);
iZc=FeatIndVec(5);
iFiltMAV=FeatIndVec(6);
iFiltMedFreq=FeatIndVec(7);
iFiltMeanFreq=FeatIndVec(8);
iFiltSsc=FeatIndVec(9);
iFiltZc=FeatIndVec(10);

MovAvgTime=0.5;
MovAvgLength=round(MovAvgTime*StimFreq);
IndexVec2=[6 8 9 10 11];
IndNum=length(IndexVec2);
IndLabels={'Unfilt','CombvEMG','CombMwave','GsvEMG','GsMwave'};

lpFiltForce = designfilt('lowpassiir','FilterOrder',8, ...
    'PassbandFrequency',5,'PassbandRipple',0.001, ...
    'SampleRate',1/S.Ts);

lpFiltMAV = designfilt('lowpassiir','FilterOrder',6, ...
    'PassbandFrequency',2,'PassbandRipple',0.01, ...
    'SampleRate',StimFreq);

for iExp=1:ExpNum
    ExpLabel=P.Labels.ExpLabels{iExp};
    
    TrialNum=P.(ExpLabel).TrialNum;
    
    for iTrial=1:TrialNum
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        x1=P.(ExpLabel).(TrialLabel).data(:,iForce);
        x1Frames=P.(ExpLabel).(TrialLabel).ForceFrames;
        
        AvgForceFrame=movmean(mean(x1Frames'),[MovAvgLength-1 0]);
        P.(ExpLabel).(TrialLabel).AvgForceFrame=AvgForceFrame';
        
        FiltMeanForce = filtfilt(lpFiltForce,mean(x1Frames'));
        P.(ExpLabel).(TrialLabel).FiltMeanForce=FiltMeanForce';
        
        for iFilt=1:FiltNum
            FiltLabel=P.Labels.FiltLabels{iFilt};

            x1=P.(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames;
            x2=P.(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames;
            [~,FrameNum]=size(x1);

            MAV1=mean(abs(x1)); 
            MAV2=mean(abs(x2)); 

            MedFreq1=zeros(1,FrameNum);
            MeanFreq1=zeros(1,FrameNum);
            Ssc1=zeros(1,FrameNum);
            Zc1=zeros(1,FrameNum);
            MedFreq2=zeros(1,FrameNum);
            MeanFreq2=zeros(1,FrameNum);
            Ssc2=zeros(1,FrameNum);
            Zc2=zeros(1,FrameNum);

            for iFrame=1:FrameNum
                
                [MedFreq1(iFrame), MeanFreq1(iFrame)]=MedMeanFreq(x1(:,iFrame),fs);
                [MedFreq2(iFrame), MeanFreq2(iFrame)]=MedMeanFreq(x2(:,iFrame),fs);

                Ssc1(iFrame)=NumSsc(x1(:,iFrame));
                Ssc2(iFrame)=NumSsc(x2(:,iFrame));

                Zc1(iFrame)=NumZc(x1(:,iFrame));
                Zc2(iFrame)=NumZc(x2(:,iFrame));

            end

            MAV1(isnan(MAV1))=0;
            MedFreq1(isnan(MedFreq1))=0;
            MeanFreq1(isnan(MeanFreq1))=0;
            Ssc1(isnan(Ssc1))=0;
            Zc1(isnan(Zc1))=0;
            MAV2(isnan(MAV2))=0;
            MedFreq2(isnan(MedFreq2))=0;
            MeanFreq2(isnan(MeanFreq2))=0;
            Ssc2(isnan(Ssc2))=0;
            Zc2(isnan(Zc2))=0;

            FiltMAV1 = filtfilt(lpFiltMAV,MAV1);
            FiltMedFreq1 = filtfilt(lpFiltMAV,MedFreq1);
            FiltMeanFreq1 = filtfilt(lpFiltMAV,MeanFreq1);
            FiltSsc1 = filtfilt(lpFiltMAV,Ssc1);
            FiltZc1 = filtfilt(lpFiltMAV,Zc1);
            
            FiltMAV2 = filtfilt(lpFiltMAV,MAV2);
            FiltMedFreq2 = filtfilt(lpFiltMAV,MedFreq2);
            FiltMeanFreq2 = filtfilt(lpFiltMAV,MeanFreq2);
            FiltSsc2 = filtfilt(lpFiltMAV,Ssc2);
            FiltZc2 = filtfilt(lpFiltMAV,Zc2);


%                 AvgMAV=movmean(MAV,[MovAvgLength-1 0]);
%                 AvgMedFreq=movmean(MedFreq,[MovAvgLength-1 0]);
%                 AvgMeanFreq=movmean(MeanFreq,[MovAvgLength-1 0]);
%                 AvgSsc=movmean(Ssc,[MovAvgLength-1 0]);
%                 AvgZc=movmean(Zc,[MovAvgLength-1 0]);


            P.(ExpLabel).(TrialLabel).(FiltLabel).MAV(1:2,:)=[MAV1; MAV2]; % 1 vEMG 2 MWaves
            P.(ExpLabel).(TrialLabel).(FiltLabel).MedFreq(1:2,:)=[MedFreq1; MedFreq2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).MeanFreq(1:2,:)=[MeanFreq1; MeanFreq2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).Ssc(1:2,:)=[Ssc1; Ssc2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).Zc(1:2,:)=[Zc1; Zc2]; 
            
            P.(ExpLabel).(TrialLabel).(FiltLabel).FiltMAV(1:2,:)=[FiltMAV1; FiltMAV2]; % 1 vEMG 2 MWaves
            P.(ExpLabel).(TrialLabel).(FiltLabel).FiltFMedFreq(1:2,:)=[FiltMedFreq1; FiltMedFreq2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).FiltMeanFreq(1:2,:)=[FiltMeanFreq1; FiltMeanFreq2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).FiltSsc(1:2,:)=[FiltSsc1; FiltSsc2]; 
            P.(ExpLabel).(TrialLabel).(FiltLabel).FiltZc(1:2,:)=[FiltZc1; FiltZc2];


%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iMedFreq)=MedFreq; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iMeanFreq)=MeanFreq;
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iSsc)=Ssc; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iZc)=Zc; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltMAV)=FiltMAV; % each column a feature
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltMedFreq)=FiltMedFreq; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltMeanFreq)=FiltMeanFreq;
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltSsc)=FiltSsc; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltZc)=FiltZc; 
%             S.(ExpLabel).(TrialLabel).(FiltLabel).FeatLabels=FeatLabels;
%             S.(ExpLabel).(TrialLabel).(FiltLabel).FeatIndVec=FeatIndVec;

            
        end                 
    end
end
% DroppedFrames Features

iExp=3;
ExpLabel=P.Labels.ExpLabels{iExp};
TrialNum=P.(ExpLabel).TrialNum;

for iTrial=1:TrialNum
    TrialLabel=sprintf('Trial_%d',iTrial);

    for iFilt=1:FiltNum
        FiltLabel=P.Labels.FiltLabels{iFilt};

        x1=P.(ExpLabel).(TrialLabel).DroppedFramesEMG;
        [~,FrameNum]=size(x1);

        MAV1=mean(abs(x1)); 

        MedFreq1=zeros(1,FrameNum);
        MeanFreq1=zeros(1,FrameNum);
        Ssc1=zeros(1,FrameNum);
        Zc1=zeros(1,FrameNum);


        for iFrame=1:FrameNum

            [MedFreq1(iFrame), MeanFreq1(iFrame)]=MedMeanFreq(x1(:,iFrame),fs);
            Ssc1(iFrame)=NumSsc(x1(:,iFrame));
            Zc1(iFrame)=NumZc(x1(:,iFrame));
            
        end

        MAV1(isnan(MAV1))=0;
        MedFreq1(isnan(MedFreq1))=0;
        MeanFreq1(isnan(MeanFreq1))=0;
        Ssc1(isnan(Ssc1))=0;
        Zc1(isnan(Zc1))=0;

        P.(ExpLabel).(TrialLabel).DroppedFramesFeatMat(1:5,:)=[MAV1; MedFreq1; MeanFreq1; Ssc1; Zc1]; % 1 vEMG 2 MWaves
              
    end                 
end


%% plotting EMG features

% FeatLabels={'Mav','MedFreq','MeanFreq','Ssc','Zc'};
PlotExp=[ 3 3] ; % 1-'MVCTrials' 2-'RCCurveTrials', 3-'ExpTrials', 4-'FatigueTrials'
PlotTrial=[4 4];
PlotFrame=[8*StimFreq 11*StimFreq];
PlotFeat=[ 1 1];
FrameInd=[PlotFrame(1):PlotFrame(2)];

FeatNum=length(FeatLabels);
FeatIndVec=[1 2 3];
iMAV=1;
iMedFreq=2;
iMeanFreq=3;
IndexVec2=[5 7 8 9 10];
IndNum=length(IndexVec2);
MatLabels={'Unfilt','CombvEMG','CombMwave','GsvEMG','GsMwaves'};
FeatTitleLabels={'MAV','Median Freq.','Mean Freq.','SSC','ZC'};

FeatLabels=P.Labels.FeatLabels;
ColorVec={'--k','--r'};
AvgColorVec={'k','r'};
FixedFrameLength=80;

for iExp=PlotExp(1):PlotExp(2)
    ExpLabel=ExpLabels{iExp};
    
    TrialNum=S.(ExpLabel).iTrial-1;
    
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        
        for iFeat=PlotFeat(1):PlotFeat(2)
            FeatLabel=FeatLabels{iFeat};
            FeatFiltLabel=sprintf('Filt%s',FeatLabel);
            
            for iFilt=1:2
                FiltLabel=FiltLabels{iFilt};
%                 vEMGLabel=sprintf('%svEMGFeatMat',FiltLabel);
%                 mWaveLabel=sprintf('%sMwaveFeatMat',FiltLabel);

                iLow=FrameInd(1)*FrameLength;
                iHigh=FrameInd(2)*FrameLength;

                vEMGFeat=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,PlotFrame(1):PlotFrame(2));
                mWaveFeat=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(2,PlotFrame(1):PlotFrame(2));
                
                AvgvEMGFeat=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatFiltLabel)(1,PlotFrame(1):PlotFrame(2));
                AvgmWaveFeat=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatFiltLabel)(2,PlotFrame(1):PlotFrame(2));

                TimeFrames=P.(ExpLabel).(TrialLabel).TimeFrames(1:FixedFrameLength,PlotFrame(1):PlotFrame(2));
                Time=mean(TimeFrames);
                
%                 Time=S.(ExpLabel).(TrialLabel).data(iLow:iHigh,iTime);
                figure(iFeat)
                subplot(2,1,1)
                plot(Time,vEMGFeat,ColorVec{iFilt},'LineWidth',1)
                hold on
                plot(Time,AvgvEMGFeat,AvgColorVec{iFilt},'LineWidth',2)
                ttl=sprintf('vEMG %s at %s, TrialNum: %d, Time(s): %.2f-%.2f',FeatLabel,ExpLabel,iTrial,PlotFrame(1)/StimFreq,PlotFrame(2)/StimFreq);
                title(ttl);
                xlabel('Frames')
                ylabel(FeatLabel)
                legend({FiltLabels{1} FiltLabels{1} FiltLabels{2} FiltLabels{2}})

                subplot(2,1,2)
                plot(Time,mWaveFeat,ColorVec{iFilt},'LineWidth',1)
                hold on
                plot(Time,AvgmWaveFeat,AvgColorVec{iFilt},'LineWidth',2)

                ttl=sprintf('mWave %s at %s, TrialNum: %d, Time(s): %.2f-%.2f',FeatLabel,ExpLabel,iTrial,PlotFrame(1)/StimFreq,PlotFrame(2)/StimFreq);
                title(ttl);
                xlabel('Frames')
                ylabel(FeatLabel)
                legend({FiltLabels{1} FiltLabels{1} FiltLabels{2} FiltLabels{2}})

            end
        end                 
    end
end
%% Plot Dropped vs Filtered
clear IndTrials vMVC sMVC IndS IndV DroppedFramesIndexFixed   DroppedFrames UnfiltFeat FiltMAV vEMGFeat DroppedFeat DroppedFramesIndexFixed DroppedFrameInd x
close all 
iExp=3;
ExpLabel=P.Labels.ExpLabels{iExp};
PlotVoli=4;
PlotStim=3;
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 40];
[VoliMVCTrials]=P.(ExpLabel).RepTableMat(:,3); % 3-VoliMVC % 4-StimMVC
[StimMVCTrials]=P.(ExpLabel).RepTableMat(:,4); % 3-VoliMVC % 4-StimMVC
vMVC=VoliMVCLevels(PlotVoli);
sMVC=StimMVCLevels(PlotStim);
IndV=find(VoliMVCTrials==vMVC);
IndS=find(StimMVCTrials==sMVC);

for iV=1:length(IndV)
    x(IndS==IndV(iV))=1;
    IndTrials=IndS(x==1);
end

MarginFromDropped=5;  % frames
PlotRange1=[ 5 11]; 
PlotRange2=[ 15 19]; 
FeatInd=[ 1 4 5];

PlotRangeFrames1=PlotRange1*StimFreq;  
PlotRangeFrames2=PlotRange2*StimFreq;        
PlotRangeInd=[PlotRangeFrames1(1):PlotRangeFrames1(2),PlotRangeFrames2(1):PlotRangeFrames2(2)];
FeatNum=length(FeatInd);
PreDroppedRange= [10.8 11 ];
PostDroppedRange=[11.2 11.5]; %% based on the dropped indices 
DroppedRange1=[5.1 10.8]*StimFreq;
DroppedRange2=[16.1 18.9]*StimFreq;
boolean DroppedFrameInd;

clear FiltRelatedFramesMean FiltRelatedFramesSD UnfiltRelatedFramesMean UnfiltRelatedFramesSD  PercentDrops
                clear UnfiltRelatedFrames FiltRelatedFrames

for iFeat=1:1
    FeatLabel=P.Labels.FeatLabels{FeatInd(iFeat)};

    for iFilt=1:FiltNum-1
        FiltLabel=P.Labels.FiltLabels{iFilt};

        for iPlot=1:length(IndTrials)
            TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));

            DroppedFramesIndexFixed=P.(ExpLabel).(TrialLabel).DroppedFramesIndexFixed;
            DroppedFrameInd1=PlotRangeFrames1(1)<= DroppedFramesIndexFixed & DroppedFramesIndexFixed<=PlotRangeFrames1(2);
            DroppedFrameInd2=PlotRangeFrames2(1)<= DroppedFramesIndexFixed & DroppedFramesIndexFixed<=PlotRangeFrames2(2);
            DroppedFrameInd=DroppedFrameInd1+DroppedFrameInd2;
            DroppedFrameInd=logical(DroppedFrameInd);

            vEMGFeat(iPlot,:)=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,PlotRangeInd);
            DroppedFrames(iPlot,:)=P.(ExpLabel).(TrialLabel).DroppedFramesFeatMat(FeatInd(iFeat),DroppedFrameInd);
            UnfiltFeat(iPlot,:)=P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,PlotRangeInd);
%             FiltMAV(iPlot,:)=P.(ExpLabel).(TrialLabel).(FiltLabel).FiltMAV(1,:); % 1 vEMG 2 MWaves
            figure(1)
            subplot(2,1,iFilt)
            for iDropped=1:length(DroppedFramesIndexFixed(DroppedFrameInd))
                DroppedLabel=sprintf('Dropped_%d',iDropped);
                
                DroppedInd=DroppedFramesIndexFixed(iDropped);
                FiltRelatedFramesMean(iPlot,iDropped)=mean(P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                FiltRelatedFramesSD(iPlot,iDropped)=std(P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                UnfiltRelatedFramesMean(iPlot,iDropped)=mean(P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                UnfiltRelatedFramesSD(iPlot,iDropped)=std(P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped));
                
                UnfiltRelatedFrames.(DroppedLabel){iPlot}=P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);
                FiltRelatedFrames.(DroppedLabel){iPlot}=P.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);

                g1=P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped);
                g2=P.(ExpLabel).(TrialLabel).DroppedFramesFeatMat(FeatInd(iFeat),DroppedFrameInd);
                
                plot(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped,ones(length(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped),1)*UnfiltRelatedFramesMean(iPlot,iDropped),...
                    'LineWidth',3,'Color','k')
                hold on
            end
%           
            PercentDrops=(UnfiltRelatedFramesMean-DroppedFrames)./UnfiltRelatedFramesMean*100;

            plot(PlotRangeInd,vEMGFeat(iPlot,:),'.', 'Color','r','LineWidth',0.01)
            plot(PlotRangeInd,UnfiltFeat(iPlot,:),'.','Color','k','LineWidth',0.01)
            plot(DroppedFramesIndexFixed(DroppedFrameInd),DroppedFrames(iPlot,:),'o','Color','b','LineWidth',1.5)
            TrialString=num2str(IndTrials');
            ttl = sprintf('%s filtered vEMG  %s vs Unfiltered at %s, Stim: %d%%, Voli: %d%%, Trials: %s',FiltLabel, FeatLabel,ExpLabel,sMVC,vMVC,TrialString);
            title(ttl);
            legend('Dropped Frame Neighbourhood','Location','northwest')
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'k', 'DisplayName', sprintf('Unfiltered MAV'))
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'r', 'DisplayName', sprintf('%s Filtered MAV',FiltLabel))
            plot([NaN NaN], [NaN NaN],'o', 'Color', 'b', 'DisplayName', sprintf('Dropped Frames'))
            ylabel('EMG MAV');
            xlabel('Frames');
            ylim([0 4*10^-4])
        end

        figure
        gLabels={'Unfilt' 'Filtered' 'Dropped'};
        for iDropped=1:length(DroppedFramesIndexFixed(DroppedFrameInd))
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
                y3(iPlot)=DroppedFrames(iPlot,iDropped);
                g3=cell(1,length(y3));
                g3(:)=gLabels(3);
                g=[ g1 g2 g3];
                y=[y1 y2 y3];

            end
            [p(iDropped),t,stats]= anova1(y,g,'off');
            subplot(1,length(DroppedFramesIndexFixed(DroppedFrameInd)),iDropped)
            results{iDropped}=multcompare(stats);
            title([])
            ylabel([])
            xlabel([])

        end


% for iUnfilt=1:length(P.(ExpLabel).(TrialLabel).Unfilt.(FeatLabel)(1,-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped))
% 
%     GROUP{iDropped2+iUnfilt}=sprintf('unfilt') ;
% end

% for iDropped=1:length(DroppedFramesIndexFixed(DroppedFrameInd))
%     
%     y=[UnfiltRelatedFramesMean(:,iDropped); DroppedFrames(:,iDropped) ;FiltRelatedFramesMean(:,iDropped)];
%     g={'Unfilt' 'Unfilt' 'Unfilt' 'Unfilt' 'Unfilt' 'Unfilt' ...
%         'Dropped' 'Dropped' 'Dropped' 'Dropped' 'Dropped' 'Dropped'...
%         'Filtered' 'Filtered' 'Filtered' 'Filtered' 'Filtered' 'Filtered'};
%     [p,t,stats]= anova1(y,g,'off');
%     figure(iDropped)
%     multcompare(stats);
% end
% 

%Plotting

        for iPlot=1:length(IndTrials)
            TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));
            
            MeanDrop=mean( PercentDrops);
            StdDrop=std( PercentDrops);

            figure(10)
            subplot(2,1,iFilt)
            for iDropped=1:length(DroppedFramesIndexFixed(DroppedFrameInd))
                DroppedInd=DroppedFramesIndexFixed(iDropped);

                plot(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped,ones(length(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped),1)*UnfiltRelatedFramesMean(iPlot,iDropped),...
                'LineWidth',2,'Color','k')
                hold on
                plot(DroppedFramesIndexFixed(DroppedFrameInd),DroppedFrames(iPlot,:),'o','Color','b','LineWidth',1.5)

    %             text(DroppedInd-10,2.2*max(UnfiltRelatedFramesMean(1,:)),sprintf('M:%.2f%%',MeanDrop(iDropped)))
    %             text(DroppedInd-10,2*max(UnfiltRelatedFramesMean(1,:)),sprintf('SD:%.2f%%',StdDrop(iDropped)))
                p_vals=results{iDropped};
                text(DroppedInd-10,2.2*max(UnfiltRelatedFramesMean(1,:)),sprintf('p1: %.2f',p_vals(1,6)))
                text(DroppedInd-10,2*max(UnfiltRelatedFramesMean(1,:)),sprintf('p2: %.2f',p_vals(2,6)))
                text(DroppedInd-10,1.8*max(UnfiltRelatedFramesMean(1,:)),sprintf('p3: %.2f',p_vals(3,6)))

                ttl = sprintf('P-values(95%% CL) at %s, Stim: %d%%, Voli: %d%%, Trials: %s',ExpLabel,sMVC,vMVC,TrialString);
                title(ttl);
                legend('Dropped Frame Neighbourhood','Location','west')
                plot([NaN NaN], [NaN NaN],'o', 'Color', 'b', 'DisplayName', sprintf('Dropped Frames'))
                ylabel('EMG MAV');
                xlabel('Frames');
                ylim([0 2.6*max(UnfiltRelatedFramesMean(1,:))])

            end
        end
        text(120,2.2*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(1,1)},gLabels{p_vals(1,2)}))
        text(120,2.0*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(2,1)},gLabels{p_vals(2,2)}))
        text(120,1.8*max(UnfiltRelatedFramesMean(1,:)),sprintf('%s & %s',gLabels{p_vals(3,1)},gLabels{p_vals(3,2)}))
    end
end

%% Fatigue Analysis
% Frequency Shift
ExpLabel='FatigueTrials';
CondLabels={'Volitional','Stimulated','StimVoli'};
FeatLabels={'MAV','MedFreq','MeanFreq','Ssc','Zc'};
MatLabels={'Unfilt','CombvEMG','CombMwave','GsvEMG','GsMwaves'};
FeatNum=length(FeatLabels);
FeatInd=[1 2 3 4 5];
TrialNum=S.(ExpLabel).iTrial-1;
iVoli=1;
iStim=2;
iStimVoli=3;

% VoliStart=7;
% VoliEnd=15;
% StimStart=22;
% StimEnd=30;
% StimVoliStart=37;
% StimVoliEnd=45;
% IndVoli=[VoliStart*StimFreq:VoliEnd*StimFreq];
% IndStim=[StimStart*StimFreq:StimEnd*StimFreq];
% IndStimVoli=[StimVoliStart*StimFreq:StimVoliEnd*StimFreq];
% Ind=[IndVoli;IndStim;IndStimVoli];


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

%% plotting Unfiltered

for iFeat=1:FeatNum 
    FeatLabel2=FeatLabels{iFeat};
    
    FeatLabel='UnfiltFeatMat';

    AvgTrials=S.(ExpLabel).(FiltLabel).(FeatLabel2).AvgTrials;
    PercentShift=S.(ExpLabel).(FiltLabel).(FeatLabel2).PercentShift;
    Mdlr_sqr=S.(ExpLabel).(FiltLabel).(FeatLabel2).Mdlr_sqr;
    p_val=S.(ExpLabel).(FiltLabel).(FeatLabel2).p_val;

    figure(2)
    subplot(1,FeatNum,iFeat)
    plot(AvgTrials','LineWidth',2)
    xlbl=sprintf('Unfiltered %s',FeatLabel2);
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

figure(2)
subplot(1,FeatNum,1)
legend(CondLabels)


%% Plot Force During Fatigue




% Linear model


%% PLotting ForceVsFeat Correlation during 
iVoli=1;
iStim=2;
iStimVoli=3;
VoliStart=7;
VoliEnd=15;
StimStart=22;
StimEnd=30;
StimVoliStart=37;
StimVoliEnd=45;
IndVoli=[VoliStart*StimFreq:VoliEnd*StimFreq];
IndStim=[StimStart*StimFreq:StimEnd*StimFreq];
IndStimVoli=[StimVoliStart*StimFreq:StimVoliEnd*StimFreq];
Ind=[IndVoli;IndStim;IndStimVoli];
FiltLabel='Comb'; % anyfilter


for iCond=1:3
    CondLabel=CondLabels{iCond};
    
    for iFeat=1:FeatNum
        FeatPlotLabel=FeatPlotLabels{iFeat};
        FeatLabel=FeatLabels{iFeat};

        UnfiltFeat=S.(ExpLabel).(FiltLabel).(FeatLabel).AvgTrials(iCond,:);
        ForceCorr= S.(ExpLabel).(FiltLabel).(FeatLabel).ForceCorr(iCond);
        AvgForce= S.(ExpLabel).(FiltLabel).(FeatLabel).AvgForce(iCond,:);
        ForMdlr_sqr=S.(ExpLabel).(FiltLabel).(FeatLabel).ForMdlr_sqr(iCond);
        Forp_val=S.(ExpLabel).(FiltLabel).(FeatLabel).Forp_val(iCond);

        figure(iCond)
        subplot(1,FeatNum,iFeat)
        plot(UnfiltFeat,AvgForce,'-ok','LineWidth',2)
        ttl1=sprintf("Unfiltered %s EMG %s Trial: %d",CondLabel, FeatPlotLabel,iTrial);
    %     ttl2=sprintf('R-sqr:%.2f p-val:%.2f Corr:%.2f ',Mdlr_sqr,p_val,ForceCorr);
        ttl2=sprintf("R-sqr:%.2f p-val:%.2f Corr:%.2f",ForMdlr_sqr,Forp_val,ForceCorr);

        title({ttl1;ttl2});
        xlabel(FeatPlotLabel)
        ylabel('Force(N)')
        grid
    end          
end
%%
PlotTrial=1;
ExpLabel='FatigueTrials';
CondLabel=CondLabels{2};

for iFeat=1:FeatNum
    FeatPlotLabel=FeatPlotLabels{iFeat};
    FeatLabel2='UnfiltFeatMat';
    FeatLabel=FeatLabels{iFeat};

    
    for iTrial=PlotTrial:PlotTrial
        TrialLabel=sprintf('Trial_%d',iTrial);

        avg_x=S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatLabel2)(Ind(iCond,:),iFeat+FeatNum);
        avg_y=S.(ExpLabel).(TrialLabel).AvgForceFrame(Ind(iCond,:),:);

        Mdlr_sqr=S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel).(CondLabel).Mdlr_sqr;
        Mdlcoef=S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel).(CondLabel).Mdlcoef;
        p_val=S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel).(CondLabel).p_val;
        iCond=S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel).(CondLabel).iCond;
        FeatForceCorr=S.(ExpLabel).(TrialLabel).AvgForceMdl.(FeatLabel).(CondLabel).FeatForceCorr;

        figure(2)
        subplot(1,FeatNum,iFeat)
        plot(avg_x,avg_y,'k','LineWidth',2)
        ttl1=sprintf('Unfiltered Stim only EMG %s Trial: %d', FeatPlotLabel,iTrial);
        ttl2=sprintf('R-sqr:%.2f p-val:%.2f Corr:%.2f',Mdlr_sqr,p_val,FeatForceCorr);
        title({ttl1;ttl2});
        xlabel(sprintf('Moving Avg (200ms) %s',FeatPlotLabel))
        ylabel('Force(N)')

    end
end          

%% Occlusion analysis
% Smooting filters
ExpLabel='ExpTrials';
MovAvgTime=0.015;
MovAvgLength=MovAvgTime*fs;
TagetLine=[0 1];
iMovAvgForce=11;
iFiltFiltForce=12;
iFiltFiltMAV=13;

lpFiltForce = designfilt('lowpassiir','FilterOrder',8, ...
    'PassbandFrequency',5,'PassbandRipple',0.001, ...
    'SampleRate',1/S.Ts);

TrialNum=S.(ExpLabel).iTrial-1;
for iTrial=1:TrialNum
    TrialLabel=sprintf('Trial_%d',iTrial);

    Force=S.(ExpLabel).(TrialLabel).data(:,iForce);

    y1 = filtfilt(lpFiltForce,Force);
    y2=movmean(Force,[MovAvgLength-1 0 ]);

    S.(ExpLabel).(TrialLabel).data(:,iFiltFiltForce)=y1;
    S.(ExpLabel).(TrialLabel).data(:,iMovAvgForce)=y2;


end



%%

subplot(2,1,1)
plot([y1 y2 MovAvgForce])
title("Filtered Waveforms")
legend("Zero-phase Filtering","Conventional Filtering")
StimMVCLevels=[ 0 10 20 30 ];
VoliMVCLevels=[10 20 30 40];


%% 
% Average reduction during stim-off

clear x
ExpLabel='ExpTrials';
PlotVoli=1;
PlotStim=4;
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 40];
ExpMat=S.(ExpLabel).RepTableMat(:,3:5); % 1-VoliMVC % 2-StimMVC 3-PW
VoliMVCTrials=ExpMat(:,1);
StimMVCTrials=ExpMat(:,2);
vMVC=VoliMVCLevels(PlotVoli);
sMVC=StimMVCLevels(PlotStim);
IndV=find(VoliMVCTrials==vMVC);
IndS=find(StimMVCTrials==sMVC);

for iV=1:length(IndV)
    x(IndS==IndV(iV))=1;
    IndTrials=IndS(x==1);
end

PlotRange=[ 9 12];
PreStimOffRange= [10.8 11 ];
PostStimOffRange=[11.2 11.3];
StimOff=11;

close all
clear PostStimOffForce PreStimOffForce
for iPlot=1:length(IndTrials)
    
    TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));
    FiltForce=S.(ExpLabel).(TrialLabel).data(:,iFiltFiltForce);
    Time=S.(ExpLabel).(TrialLabel).data(:,iTime);
    PreStimOffInd=PreStimOffRange(1)<=Time & Time<=PreStimOffRange(2);
    PostStimOffInd=PostStimOffRange(1)<=Time & Time<=PostStimOffRange(2);

    figure(1)
    subplot(3,1,1)
    plot(Time,FiltForce)
    hold on
    lgd{iPlot}=sprintf('Trial %d',IndTrials(iPlot));
    
    figure(2)
    PlotInd=PlotRange(1)<=Time & Time<=PlotRange(2);
    plot(Time(PlotInd),FiltForce(PlotInd))
    hold on
    lgd{iPlot}=sprintf('Trial %d',IndTrials(iPlot));
    
    PreStimOffForce(iPlot)=mean(FiltForce(PreStimOffInd ));
    PostStimOffForce(iPlot)=mean(FiltForce(PostStimOffInd ));

end


MeanPreStimOffForce=mean(PreStimOffForce);
MeanPostStimOffForce=mean(PostStimOffForce);
PercentDrop=(PreStimOffForce-PostStimOffForce)./PreStimOffForce*100;
MeanPercentDrop=mean(PercentDrop);
SDPercentDrop=std(PercentDrop);
text(9.2,MeanPreStimOffForce*1.24,'Percent Drops:','FontSize',14)
text(9.2,MeanPreStimOffForce*1.2,sprintf('%.2f%%, ',PercentDrop),'FontSize',14)
text(9.2,MeanPreStimOffForce*1.15,sprintf('Mean:%.2f%%, SD:%.2f%%',MeanPercentDrop,SDPercentDrop),'FontSize',14)

plot([ PreStimOffRange(1) PreStimOffRange(2)], [mean(PreStimOffForce) mean(PreStimOffForce)] ,'Color','k','LineWidth',1.5);
plot([ PostStimOffRange(1) PostStimOffRange(2)], [mean(PostStimOffForce) mean(PostStimOffForce)] ,'Color','k','LineWidth',1.5);

plot([PreStimOffRange(1) PreStimOffRange(1)],[0 100],'--','Color','r')
plot([StimOff StimOff],[0 100],'--','Color','r')
plot([PostStimOffRange(1) PostStimOffRange(1)],[0 100],'--','Color','r')
plot([PostStimOffRange(2) PostStimOffRange(2)],[0 100],'--','Color','r')

ylim([MeanPreStimOffForce/2 MeanPreStimOffForce*1.3])
xlabel('Time(s)')
ylabel('Force(N)')
grid 
ttl1=sprintf('Force at %s, VoliMVC: %d%% StimMVC: %d%% Trials:',ExpLabel,vMVC,sMVC);
ttl2=sprintf(' ,%d',IndTrials);
title(strcat(ttl1,ttl2))
legend(lgd)

PW=S.(ExpLabel).(TrialLabel).data(:,iPW);
NormPW= max(FiltForce)*PW/max(PW)/2;
figure(1)
plot(Time,NormPW)
% lgd2={lgd ' Avg. Force', 'Stim Off Time'};
ttl1=sprintf('Force at %s, VoliMVC: %d%% StimMVC: %d%% Trials:',ExpLabel,vMVC,sMVC);
ttl2=sprintf(' ,%d',IndTrials);
title(strcat(ttl1,ttl2))
xlabel('Time(s)')
ylabel('Force(N)')
legend(lgd)
grid

figure(1)
vEMGLabels={'CombvEMG', 'GsvEMG','Unfilt'};
for iPlot=1:length(IndTrials)
    TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));

    for  iFilt=1:FiltNum
    FiltLabel=FiltLabels{iFilt};
    FeatMatLabel=sprintf('%sFeatMat',vEMGLabels{iFilt});
    
    FiltMAV=S.(ExpLabel).(TrialLabel).(FiltLabel).(FeatMatLabel)(:,iFiltMAV); % each column a feature
    Time=(1:length(FiltMAV))/StimFreq;
    
    figure(1)
    subplot(3,1,iFilt+1)
    plot(Time,FiltMAV)
    hold on
    lgd{iPlot}=sprintf('Trial %d',IndTrials(iPlot));
    ttl1=sprintf('%s vEMG MAV at %s, VoliMVC: %d%% StimMVC: %d%% Trials:',FiltLabel,ExpLabel,vMVC,sMVC);
    ttl2=sprintf(' ,%d',IndTrials);
    title(strcat(ttl1,ttl2))
    xlabel('Time(s)')
    ylabel(sprintf('%s estimated MAV',FiltLabel));
    grid on
    end
end



%%
function [Mwave, vEMG]=FiltComb(x)
    % performs comb filter 
    % x: EMG signal where columns are frames and rows are samples
   [~, FrameNum] =size(x);

    for iFrame=2:FrameNum
        y(:,iFrame) = (x(:,iFrame)-x(:,iFrame-1))/sqrt(2);

    end

    vEMG=y;
    Mwave=x-y;

end


function filtered_sig = SUB_GS_filter(input_sig, L, M)
%GS_FILTER Summary of this function goes here
%   Detailed explanation goes here

%L = frame length
%M = filter order

[r,c] = size(input_sig);
epsilon = zeros(r,c,M);
epsilon(:,:,1) = fliplr(input_sig);
initial_sig = epsilon(:,1,1);

for m = (0:1:M-1)+1
    for i = (0:1:M-m-1)+1
%         w(i,m) = (epsilon(:,i,m)'*epsilon(:,M-m+1,m))/(sum(epsilon(:,M-m+1,m).^2));
        w(i,m) = (epsilon(:,i,m)'*epsilon(:,M-m+1,m))/(norm(epsilon(:,M-m+1,m))^2);
        if (sum(epsilon(:,M-m+1,m).^2)) == 0
            disp('NaN detected')
        end
        epsilon(:,i,m+1) = epsilon(:,i,m) - w(i,m)*epsilon(:,M-m+1,m);
    end
end

filtered_sig = epsilon(:,:,end);
filtered_sig(:,2) = initial_sig-epsilon(:,1,end);

% if isnan(filtered_sig(1,2))
%     disp('DEBUG')
% end

end

function zc=NumZc(x)
zc=0;
thres=.000004;

for i=2:length(x)
    if x(i-1)*x(i)<0 && abs(x(i-1)-x(i))>thres
        zc=zc+1;
    end
end

end
function SSC=NumSsc(X)
% x vector 

% thres=.0000025;
thres=.0000035;
N=length(X); SSC=0;

for i=2:N-1
  if ((X(i) > X(i-1) && X(i) > X(i+1)) || (X(i) < X(i-1) && X(i) < X(i+1))) ...
      && ((abs(X(i)-X(i+1)) >= thres) || (abs(X(i)-X(i-1)) >= thres))
    SSC=SSC+1; 
  end

end

end

function [MedFreq MeanFreq]=MedMeanFreq(x,fs)
%calculates mean and med freq of a signal using psd method

xdft = fft(x);
xdft = xdft(1:floor(length(x)/2)+1);
xdft(2:end-1) = 2*xdft(2:end-1);
psd = 1/(length(x)*fs)*abs(xdft).^2;
freq = 0:fs/length(x):fs/2;


MedFreq = medfreq(psd,freq);
MeanFreq = meanfreq(psd,freq);

end