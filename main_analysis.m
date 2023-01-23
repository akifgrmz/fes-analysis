%% main_analysis v0.1
% this file is for running the analysis for one file only but will be a
% function for running all the files in future
%% Run tidy_data if you have not yet done so 

clear all
FolderName="jan12";
tidy_data(FolderName);

%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12" ];

for iTest=1:length(TestFolders)
    TestFiles{iTest}=sprintf("%s_test",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%%Defining Initial Parameters
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    fs=S.(TestStruct).ExpPar.fs;
    
    %%Defining the parameters 
% Experiment Parameters

%     StimTime=0.007;
%     TrigDelay=0.006;
%     TrigTime=0.003;

    % Experiments indices 
    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    DataInd=S.(TestStruct).ExpPar.DataInd;
    TableInd=DataInd.Properties.VariableNames;
    iForce=table2array(DataInd(:,"Force"));
    iTrigger=table2array(DataInd(:,"Trigger"));
    iEMG=table2array(DataInd(:,"EMG"));
    iTime=table2array(DataInd(:,"Time"));
    iPW=table2array(DataInd(:,"PW"));

    S.(AnaStruct).AnaPar.ExpTable=table(ExpLabels(1),ExpLabels(2),...
        ExpLabels(3),ExpLabels(4),ExpLabels(5),'VariableNames',["MVC","RC","Cus","Occ","Fat"]);

    % Analysis Parameters
    S.(AnaStruct).AnaPar.AnaLabels=["BPFilt_EMG"];
    S.(AnaStruct).AnaPar.AnaInd=table([],'VariableNames',S.(AnaStruct).AnaPar.AnaLabels);
    
    % TrigLowThres=0.5;
    % TrigHighThres=4;
    % TrigThres=1;

    DataLabels= {'EMG', 'Trigger Signal', 'Force (N)','Pulse Width','Time (s)'...
        ,'BP Filtered EMG','Blanked EMG','Comb Filtered vEMG','Comb Filtered m-Waves'...
        ,'GS Filtered vEMG','GS Filtered m-Waves'};
    
    FeatLabels={'MAV','MedFreq','MeanFreq','Ssc','Zc'};
    FiltLabels={'Unfilt','Comb','GS'};
    S.(AnaStruct).AnaPar.DataLabels=DataLabels;
    S.(AnaStruct).AnaPar.FeatLabels=FeatLabels;
    S.(AnaStruct).AnaPar.FiltLabels=FiltLabels;


    %Incorporate the redo trials at the begginning 
    %Design 10-500Hz 20th order butterworth for filtfilt

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
    S.(AnaStruct).AnaPar.BPFilter.d1=d1;
end

%%Incorporate Redo Trials 

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpstoAna=ExpLabels([2,4,5]);

    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        
        RedoTrials=S.(TestStruct).(ExpLabel).RedoTrials;
        
        if ~isempty(RedoTrials)
           NumofRedos=length(RedoTrials);
           
           for iRedo=1:NumofRedos
              
                RedoLabel=sprintf('RedoTrial_%d',RedoTrials(iRedo));
                TrialLabel=sprintf('Trial_%d',RedoTrials(iRedo));

                S.(TestStruct).(ExpLabel).(TrialLabel)=S.(TestStruct).(ExpLabel).(RedoLabel);
               
           end
        end 
    end
end

%

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    AmpGain=S.(TestStruct).ExpPar.AmpGain;

    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels{iExp};

        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
  
            x=S.(TestStruct).(ExpLabel).(TrialLabel).data.(TableInd{iEMG})/AmpGain;
          
            BPFilt_EMG = filtfilt(d1,x);
            S.(AnaStruct).(ExpLabel).(TrialLabel).data=table(BPFilt_EMG);
        end
    end
end

%Trigger and Blanking

BlankTime=0.0035;
BlankLength=round(BlankTime*fs)+1;
ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
ExpstoAna=ExpLabels([2,3,4,5]);

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    S.(AnaStruct).AnaPar.BlankLength=BlankLength;
    S.(AnaStruct).AnaPar.BlankTime=BlankTime;
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);

            y=S.(AnaStruct).(ExpLabel).(TrialLabel).data.("BPFilt_EMG");
            x=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Trigger");

            Time=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Time"); 
                % # blanking
            RisingInd=diff(x)>0;
            FallingInd=diff(x)<0;
            BegofFrames=find(FallingInd>0);
            FallingInd(end+1)=0;
            RisingInd(end+1)=0;

            FrameLength=BegofFrames(11)-BegofFrames(10);  % could be any two frames
            KeepLength=FrameLength-BlankLength;
            z=zeros(length(y),1);
            for iFall=2:length(BegofFrames)

                z(BegofFrames(iFall)-KeepLength:BegofFrames(iFall))=y(BegofFrames(iFall)-KeepLength:BegofFrames(iFall));
            end
    %         
            S.(AnaStruct).(ExpLabel).(TrialLabel).Indices=table([FallingInd], [RisingInd]);
            S.(AnaStruct).(ExpLabel).(TrialLabel).data.("BlankedEMG")=z;
            S.(AnaStruct).(ExpLabel).(TrialLabel).data.("Time")=Time;
            S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames=BegofFrames;
            S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength=FrameLength;

        end
    end
end

%% Plotting after blanking
Lbl='Occ';
PlotTrial=[7 7 ];
PlotFrame=floor([ 487 489]);

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
    ExpLabel=ExpTable{1};
    DataLabels=S.(AnaStruct).AnaPar.DataLabels;
    AmpGain=S.(TestStruct).ExpPar.AmpGain;

    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);

        FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;
        BegofFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
    %     Ind=[BlankLength+PlotFrame(1)*FrameLength BlankLength+(PlotFrame(2))*FrameLength];
        Ind= [BegofFrames(PlotFrame(1)) BegofFrames(PlotFrame(2)+1)];

        BPFilt_EMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BPFilt_EMG');
        BlankEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BlankedEMG');
        Trig=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Trigger');
        EMG=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('EMG');
        Time=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Time');

        figure(1)
        subplot(length(TestFolders),1,iTest)
        plot(Time,Trig/1000,'b','LineWidth',2)
        hold on
        plot(Time,EMG/AmpGain,'k','LineWidth',2)
        plot(Time,BlankEMG,'r','LineWidth',2)
        legend({'Trigger','Pre-Blanking and Filtering', 'Post-Blanking and Filtering'})
        ttl=sprintf('BP Filtered %s at %s,Test: %s, TrialNum: %d, Frame: %d-%d,'...
            ,DataLabels{1},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2));
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        ylim([1.5*min(BPFilt_EMG) 1.5*max(BPFilt_EMG)])

    end
end

%%Frames Matrix for Force and EMG
% (also remove blanked periods)
% (Add unfilt as a filter and add e m-wave as zero matrix)
% Identify Zero stim trials -> Done

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpstoAna=ExpLabels([2,3,4,5]);
    
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        
        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            
            PW=S.(TestStruct).(ExpLabel).(TrialLabel).data.('PW');
            BegofFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);

            clear PWofFrames
            for iFrame=1:NumofFrames
                %PW
                PWofFrames(iFrame)=PW(BegofFrames(iFrame));
            end

            S.(AnaStruct).(ExpLabel).(TrialLabel).PWofFrames=PWofFrames';
        end
    end

    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpstoAna=ExpLabels([2,3,4,5]);
    BlankLength=S.(AnaStruct).AnaPar.BlankLength;
    
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};

        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;
            
            BlankEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data.('BlankedEMG');
            BPFilt_EMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data.('BPFilt_EMG');

            Force=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Force');
            Time=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Time');

            BegofFrames= S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);

            clear y x f t FrameLengths PWFrames TmFrames
            for iFrame=1:NumofFrames-1
                FrameLengths(iFrame)=-BegofFrames(iFrame)+BegofFrames(iFrame+1);
                % # BP filt EMG
                x(:,iFrame)=BPFilt_EMG(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
                % # Blanked EMG
                y(:,iFrame)=BlankEMG(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
                % # Force
                f(:,iFrame)= Force(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
                % # Time Frames
                t(:,iFrame)= Time(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));

            end

            S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames=y;
            S.(AnaStruct).(ExpLabel).(TrialLabel).EMGFrames=x;
            S.(AnaStruct).(ExpLabel).(TrialLabel).ForceFrames=f;
            S.(AnaStruct).(ExpLabel).(TrialLabel).TmFrames=t;
            S.(AnaStruct).(ExpLabel).(TrialLabel).BlankTmFrames=t;
            S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLengths=FrameLengths;

        end
    end
end

%%Extracting Dropped Frames 
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.('Occ'));

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
    TargetProfile=S.(TestStruct).(ExpLabel).TargetProfile2;
    stim_freq=S.(TestStruct).ExpPar.FreqList(1);
    TrialsPW=S.(TestStruct).(ExpLabel).TrialsPW;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameRange=[(stim_freq+1)*TargetProfile(1,2),stim_freq*TargetProfile(1,4)];
        PW=S.(TestStruct).(ExpLabel).(TrialLabel).data.('PW');

        BegofFrames= S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
        NumofFrames=length(BegofFrames);

        if TrialsPW(iTrial) == 0
            ZeroInd=zeros(NumofFrames,1);
        else
            ZeroInd=find(( PW(BegofFrames)==0)==1);
        end

        DroppedFrames= ZeroInd(ZeroInd>=FrameRange(1) & ZeroInd<=FrameRange(2) );

        DroppedEMG= S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames(:,DroppedFrames);

        S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames=DroppedFrames;
        S.(AnaStruct).(ExpLabel).(TrialLabel).ZeroInd=ZeroInd;
        S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedEMG=DroppedEMG;

    end
end


%%Target Lines by samples and frames
% this can go to tidyscript

% ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
% ExpstoAna=ExpLabels([4]);
% fs= S.(TestStruct).ExpPar.fs;
% TargetLabel={"TargetProfile2"}
% for iExp=1:length(ExpstoAna)
%     ExpLabel=ExpstoAna{iExp};
%     NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
%     
%     TargetPro=S.(TestStruct).(ExpLabel).TargetProfile2
% 
%     
% %     for iTrial=1:NumofTrials
% %         TrialLabel=sprintf('Trial_%d',iTrial);
% %         
% %         y=S.(AnaStruct).(ExpLabel).(TrialLabel).data.("BPFilt_EMG");
% %         x=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Trigger");
% %         
% %         Time=S.(TestStruct).(ExpLabel).(TrialLabel).data.("Time"); 
% % 
% % 
% %     end
% end


%% Plotting the dropped frames 

iTrial=9;
TrialLabel=sprintf('Trial_%d',iTrial);
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.('Occ'));

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});

    S.(TestStruct).(ExpLabel).TrialsPW
    DataLabels=S.(AnaStruct).AnaPar.DataLabels;
    DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
    DroppedEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedEMG;
    clear lgd

    for iFrame=1:length(DroppedFrames)
        figure(iTrial)
        subplot(length(TestFolders),1,iTest)
        plot(DroppedEMG(:,iFrame),'LineWidth',2)
        hold on
        lgd(iFrame)=sprintf("Frame %d",DroppedFrames(iFrame));
        legend(lgd,'Location','NorthWest');
        ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
        title(ttl);
        xlabel('Samples')
        ylabel('BP Filtered EMG')

    end
end

%% M-wave filtering  
% (Remove blanked periods before filtering, Do unfiltered Frames)
% (update to recreate the vector forms of the signals)

gs_order=3;

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpstoAna=ExpLabels([2,3,4,5]);
    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;

    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;
        TrialsPW=S.(TestStruct).(string(S.(AnaStruct).AnaPar.ExpTable.('Occ'))).TrialsPW;

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            clear T Time TimeFrames x_frames
            BegofFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
            FrameNum=length(BegofFrames);
            
            if string(ExpLabel)==string(S.(AnaStruct).AnaPar.ExpTable.('Occ')) && TrialsPW(iTrial) ~= 0
                DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                KeepInd=setdiff(1:FrameNum-1,DroppedFrames);

                x_frames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames(:,KeepInd); % excluding dropped frames 
                [FrameLength, FrameNum]=size(x_frames);
                TimeFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankTmFrames(:,KeepInd);
                Time=reshape(TimeFrames,[1],FrameLength*FrameNum);
                S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength=FrameLength;

            else 
                x_frames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames;
                [FrameLength, FrameNum]=size(x_frames);

                KeepInd=[];
                TimeFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).TmFrames;
                Time=reshape(TimeFrames,[1],FrameLength*FrameNum);
            end


            S.(AnaStruct).(ExpLabel).(TrialLabel).KeepInd=KeepInd;
            S.(AnaStruct).(ExpLabel).Time=table(Time','VariableName',"Time");

            % # Unfilt
            iFilt=1;
            FiltLabel=FiltLabels{iFilt};
            UnfiltvEMGFrames=x_frames;
            UnfiltMWaveFrames=zeros(size(UnfiltvEMGFrames));
            UnfiltvEMG=reshape(UnfiltvEMGFrames,1,FrameLength*FrameNum);
            UnfiltMWave=reshape(UnfiltMWaveFrames,1,FrameLength*FrameNum);
            T=[table(UnfiltvEMG','VariableName',"UnfiltvEMG") table(UnfiltMWave','VariableName',"UnfiltMWave")];

            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=UnfiltMWaveFrames;
            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=UnfiltvEMGFrames;
            AllLabel=sprintf("%s_all",FiltLabel);
            S.(AnaStruct).(ExpLabel).(FiltLabel).T=T;   


            if string(ExpLabel)==string(S.(AnaStruct).AnaPar.ExpTable.('Occ')) && TrialsPW(iTrial) ~= 0

                BlankEMGFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=UnfiltvEMGFrames;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGwithDropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=UnfiltMWaveFrames;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).mWaveswithDropped=BlankEMGFrames;
            end

            % # Comb Filter
            iFilt=2;

            [CombMWaveFrames,CombvEMGFrames]=FiltComb(x_frames);
            CombvEMG=reshape(CombvEMGFrames,1,FrameLength*FrameNum);
            CombMWave=reshape(CombMWaveFrames,1,FrameLength*FrameNum);

            FiltLabel=FiltLabels{iFilt};
            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=CombMWaveFrames;
            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=CombvEMGFrames;
            T=[table(CombvEMG','VariableName',"CombvEMG") table(CombMWave','VariableName',"CombMWave")];
            AllLabel=sprintf("%s_all",FiltLabel);
            S.(AnaStruct).(ExpLabel).(FiltLabel).T=T;   

            if string(ExpLabel)==string(S.(AnaStruct).AnaPar.ExpTable.('Occ')) && TrialsPW(iTrial) ~= 0

                BlankEMGFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=CombvEMGFrames;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGwithDropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=CombMWaveFrames;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).mWaveswithDropped=BlankEMGFrames;
            end

            % # GS Filter
            % This method is taking too long to calculate
            iFilt=3;
            FiltLabel=FiltLabels{iFilt};

            clear vEMGMat MWaveMat
            for iFrame = 1:FrameNum-gs_order
                temp_vect = SUB_GS_filter(x_frames(:,iFrame:gs_order+iFrame),FrameLength,gs_order);
                vEMGMat(:,iFrame+gs_order) = temp_vect(:,1);
                MWaveMat(:,iFrame+gs_order) = temp_vect(:,2);
            end

            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=vEMGMat;
            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=MWaveMat;

            GSvEMG=reshape(vEMGMat,1,FrameLength*FrameNum);
            GSMWave=reshape(MWaveMat,1,FrameLength*FrameNum);

            T=[table(GSvEMG','VariableName',"GSvEMG") table(GSMWave','VariableName',"GSMWave")];
            AllLabel=sprintf("%s_all",FiltLabel);
            S.(AnaStruct).(ExpLabel).(FiltLabel).T=T;   

            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMG=GSvEMG;
            S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWave=GSMWave;

            if string(ExpLabel)==string(S.(AnaStruct).AnaPar.ExpTable.('Occ')) && TrialsPW(iTrial) ~= 0

                BlankEMGFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=vEMGMat;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGwithDropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=MWaveMat;
                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).mWaveswithDropped=BlankEMGFrames;
            end
            % # Blanking filter comes here

        end
    end
end


%% Plotting, Debugging 
close all
Lbl='RC';
PlotTrial=[ 12 12];
PlotTime=[6.5 6.8];
PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(Lbl));

    % PlotFrame=floor([ 300 305]); % first frame is skipped

    ExpLabels=S.(TestStruct).ExpPar.ExpLabels;
    ExpstoAna=ExpLabels([2,3,4,5]);
    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;
    clear EMGFrames MWaveFrames vEMG MWave

    NumofTrials=S.(TestStruct).(ExpLabel).NumofTrials;

    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;

        for iFilt=1:length(FiltLabels)
            FiltLabel=FiltLabels{iFilt};

            EMGFrames(:,:,iFilt)=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames(:,:,iFilt)=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
            TimeFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
            Time=reshape(TimeFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            Trigger=S.(TestStruct).(ExpLabel).(TrialLabel).data.('Trigger');
            vEMG(:,iFilt)=reshape(EMGFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave(:,iFilt)=reshape(MWaveFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

        end

        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{2},TestFolders{iTest}),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,vEMG(:,2),'r','LineWidth',2)
        hold
        legend({'Unfiltered', 'Comb vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d'...
            ,DataLabels{8},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,2),'r','LineWidth',2)
        legend({'Unfiltered', 'Comb MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{9},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
    
        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{3},TestFolders{iTest}),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,vEMG(:,3),'r','LineWidth',2)
        legend({'Unfiltered', 'GS vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{10},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})
        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,3),'r','LineWidth',2)
        legend({'Unfiltered', 'GS MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{11},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{iTime})
        ylabel(DataLabels{iEMG})

    end
end


%% EMG Features
% (move force averaging to the top where BP is held)


ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(Lbl));

% lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
%     'PassbandFrequency',1,'PassbandRipple',0.1,'SampleRate',stim_freq);

Fpass = 2;
Fstop = 10;
Ap = 0.1;
Ast = 30;
Fs = 2000; 

lpFilt = designfilt('lowpassfir','PassbandFrequency',Fpass,...
  'StopbandFrequency',Fstop,'PassbandRipple',Ap,...
  'DesignMethod', 'kaiserwin','SampleRate',stim_freq);

% lpFilt = designfilt('lowpassiir','PassbandFrequency',Fpass,...
%   'StopbandFrequency',Fstop,'PassbandRipple',Ap,'SampleRate',stim_freq);

% fvtool(lpFilt,'MagnitudeDisplay','magnitude')
N = filtord(lpFilt)
%
clear FiltFramesInd
stim_freq=S.(TestStruct).ExpPar.FreqList(1);
fs=S.(TestStruct).ExpPar.fs;
ExptoAnaInd=[2,4,5];
FiltFramesInd{ExptoAnaInd(1)}=floor([5.5*stim_freq: 10*stim_freq]);
FiltFramesInd{ExptoAnaInd(2)}=floor([5.5*stim_freq: 15*stim_freq]);
FiltFramesInd{ExptoAnaInd(3)}=floor([5.5*stim_freq: 35*stim_freq]);

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    ExptoAna=S.(TestStruct).ExpPar.ExpLabels([2,4,5]);
    for iExp=1:length(ExptoAna)
        ExpLabel=ExptoAna{iExp};
        
        TrialNum=S.(TestStruct).(ExpLabel).NumofTrials;

        for iTrial=1:TrialNum
            TrialLabel=sprintf('Trial_%d',iTrial);

            x1Frames=S.(AnaStruct).(ExpLabel).(TrialLabel).ForceFrames;

            FiltForceFrames = filtfilt(lpFilt,mean(x1Frames)); % mean force for each frame so fs=35hz
            S.(AnaStruct).(ExpLabel).(TrialLabel).FiltForceFrames=FiltForceFrames;

            for iFilt=1:length(S.(AnaStruct).AnaPar.FiltLabels)
                FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};

                vEMGFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames;
                MWaveFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames;
                [~,FrameNum]=size(vEMGFrames);

                MAV_vEMG=mean(abs(vEMGFrames)); 
                MAV_MWaves=mean(abs(MWaveFrames)); 

                MedFreq_vEMG=zeros(1,FrameNum);
                MeanFreq_vEMG=zeros(1,FrameNum);
                MedFreq_MWaves=zeros(1,FrameNum);
                MeanFreq_MWaves=zeros(1,FrameNum);

                Ssc1=zeros(1,FrameNum);
                Zc1=zeros(1,FrameNum);
                Ssc2=zeros(1,FrameNum);
                Zc2=zeros(1,FrameNum);

                for iFrame=1:FrameNum

                    [MedFreq_vEMG(iFrame), MeanFreq_vEMG(iFrame)]=MedMeanFreq(vEMGFrames(:,iFrame),fs);
                    [MedFreq_MWaves(iFrame), MeanFreq_MWaves(iFrame)]=MedMeanFreq(MWaveFrames(:,iFrame),fs);

                    Ssc1(iFrame)=NumSsc(vEMGFrames(:,iFrame));
                    Ssc2(iFrame)=NumSsc(MWaveFrames(:,iFrame));

                    Zc1(iFrame)=NumZc(vEMGFrames(:,iFrame));
                    Zc2(iFrame)=NumZc(MWaveFrames(:,iFrame));

                end

%                 MAV_vEMG(isnan(MAV_vEMG))=0;
                MedFreq_vEMG(isnan(MedFreq_vEMG))=0;
                MeanFreq_vEMG(isnan(MeanFreq_vEMG))=0;
                Ssc1(isnan(Ssc1))=0;
                Zc1(isnan(Zc1))=0;
                MAV_MWaves(isnan(MAV_MWaves))=0;
                MedFreq_MWaves(isnan(MedFreq_MWaves))=0;
                MeanFreq_MWaves(isnan(MeanFreq_MWaves))=0;
                Ssc2(isnan(Ssc2))=0;
                Zc2(isnan(Zc2))=0;

                FiltMAV_vEMG = filtfilt(lpFilt,MAV_vEMG);
                FiltMedFreq_vEMG = filtfilt(lpFilt,MedFreq_vEMG);
                FiltMeanFreq_vEMG = filtfilt(lpFilt,MeanFreq_vEMG);
                FiltMAV_MWaves = filtfilt(lpFilt,MAV_MWaves);
                FiltMedFreq_MWaves = filtfilt(lpFilt,MedFreq_MWaves);
                FiltMeanFreq_MWaves = filtfilt(lpFilt,MeanFreq_MWaves);
    %             
    %             FiltSsc2 = filtfilt(lpFiltMAV,Ssc2);
    %             FiltZc2 = filtfilt(lpFiltMAV,Zc2);
    %             FiltSsc1 = filtfilt(lpFiltMAV,Ssc1);
    %             FiltZc1 = filtfilt(lpFiltMAV,Zc1);

                S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats=...
                     table(MAV_vEMG',MedFreq_vEMG',MeanFreq_vEMG',MAV_MWaves',MedFreq_MWaves',MeanFreq_MWaves',...
                     'VariableNames',["MAV_vEMG" "MedFreq_vEMG" "MeanFreq_vEMG" "MAV_MWaves" "MedFreq_MWaves" "MeanFreq_MWaves"]);

                 S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats=...
                     table(FiltMAV_vEMG',FiltMedFreq_vEMG',FiltMeanFreq_vEMG',FiltMAV_MWaves',FiltMedFreq_MWaves',FiltMeanFreq_MWaves',...
                     'VariableNames',[ "Filt_MAV_vEMG" "Filt_MedFreq_vEMG" "Filt_MeanFreq_vEMG" "Filt_MAV_MWaves" "Filt_MedFreq_MWaves" "Filt_MeanFreq_MWaves"]);

    %             P.(ExpLabel).(TrialLabel).(FiltLabel).FiltSsc(1:2,:)=[FiltSsc1; FiltSsc2]; 
    %             P.(ExpLabel).(TrialLabel).(FiltLabel).FiltZc(1:2,:)=[FiltZc1; FiltZc2];

            end                 
        end
    end
end

%%DroppedFrames Features

Lbl='Occ';
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(Lbl));


for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    TrialNum=S.(TestStruct).(ExpLabel).NumofTrials;
    BlankLength=S.(AnaStruct).AnaPar.BlankLength;

    for iTrial=1:TrialNum
        TrialLabel=sprintf('Trial_%d',iTrial);

        for iFilt=1:length(S.(AnaStruct).AnaPar.FiltLabels)
            FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};

            DroppedEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedEMG;
            [~,FrameNum]=size(DroppedEMG);

            MAV_vEMG=mean(abs(DroppedEMG)); 

            MedFreq_vEMG=zeros(1,FrameNum);
            MeanFreq_vEMG=zeros(1,FrameNum);
            Ssc=zeros(1,FrameNum);
            Zc=zeros(1,FrameNum);


            for iFrame=1:FrameNum

                [MedFreq_vEMG(iFrame), MeanFreq_vEMG(iFrame)]=MedMeanFreq(DroppedEMG(:,iFrame),fs);
                Ssc(iFrame)=NumSsc(DroppedEMG(:,iFrame));
                Zc(iFrame)=NumZc(DroppedEMG(:,iFrame));

            end

            MAV_vEMG(isnan(MAV_vEMG))=0;
            MedFreq_vEMG(isnan(MedFreq_vEMG))=0;
            MeanFreq_vEMG(isnan(MeanFreq_vEMG))=0;
            Ssc(isnan(Ssc))=0;
            Zc(isnan(Zc))=0;

           S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat=...
               table(MAV_vEMG',MedFreq_vEMG',MeanFreq_vEMG',Ssc',Zc',...
               'VariableNames',["MAV_vEMG" "MedFreq_vEMG" "MeanFreq_vEMG" "Ssc" "Zc"]);
        end                 
    end
end
%% Saving analysis 
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    str=sprintf('%s/%s',TestFolders{iTest},AnaStruct);
    save(str,'-struct','S',TestStruct,AnaStruct)
end


%% plotting EMG features
clear vEMGFeat
close all

Lbl='RC';
PlotTrial=[1 1];
PlotTime=[5 10];
PlotFrame=[PlotTime(1)*stim_freq PlotTime(2)*stim_freq];
PlotFeat=[ 1 1 ]; % MAV only 
FrameInd=[PlotFrame(1):PlotFrame(2)];
Time=FrameInd/stim_freq;
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(Lbl));
FeatLabels=S.(AnaStruct).AnaPar.FeatLabels;
ColorVec={'--b','--k','--r'};
AvgColorVec={'b','k','r'};

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    TwitchPW=S.(TestStruct).(S.(AnaStruct).AnaPar.ExpTable.('RC')).PWPoints(1);
    
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;

        PW=S.(AnaStruct).(ExpLabel).(TrialLabel).PWofFrames(FrameInd);
        FiltForceFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).FiltForceFrames(FrameInd);

        for iFeat=PlotFeat(1):PlotFeat(2)
            FeatLabel=FeatLabels{iFeat};
            FeatFiltLabel=sprintf('Filt%s',FeatLabel);

            for iFilt=1:length( FiltLabels)
                FiltLabel=FiltLabels{iFilt};
                FiltvEMGLabel=sprintf('Filt_%s_vEMG',FeatLabel);
                FiltMWaveLabel=sprintf('Filt_%s_MWaves',FeatLabel);
                vEMGLabel=sprintf('%s_vEMG',FeatLabel);
                MWaveLabel=sprintf('%s_MWaves',FeatLabel);

                FiltvEMGFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats.(FiltvEMGLabel)(PlotFrame(1):PlotFrame(2));
                FiltmWaveFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats.(FiltMWaveLabel)(PlotFrame(1):PlotFrame(2));

                vEMGFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel)(PlotFrame(1):PlotFrame(2));
                mWaveFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(MWaveLabel)(PlotFrame(1):PlotFrame(2));

                figure(iTest)
                subplot(2,1,1)
                plot(Time,vEMGFeat,ColorVec{iFilt},'LineWidth',1)
                hold on
                plot(Time,FiltvEMGFeat,AvgColorVec{iFilt},'LineWidth',2)
                ttl=sprintf('vEMG %s at %s, TrialNum: %d, Time(s): %.2f-%.2f',FeatLabel,ExpLabel,iTrial,PlotFrame(1)/stim_freq,PlotFrame(2)/stim_freq);
                title(ttl);
                xlabel('Time')
                ylabel(FeatLabel)

                subplot(2,1,2)
                plot(Time,mWaveFeat,ColorVec{iFilt},'LineWidth',1)
                hold on
                plot(Time,FiltmWaveFeat,AvgColorVec{iFilt},'LineWidth',2)
                ttl=sprintf('mWave %s at %s, TrialNum: %d, Time(s): %.2f-%.2f',FeatLabel,ExpLabel,iTrial,PlotFrame(1)/stim_freq,PlotFrame(2)/stim_freq);
                title(ttl);
                xlabel('Time')
                ylabel(FeatLabel)
%                 legend({FiltLabels{1} FiltLabels{1} FiltLabels{2} FiltLabels{2} FiltLabels{3} FiltLabels{3}})
            end

            subplot(2,1,1)
            nrm2=100000;
            plot(Time,FiltForceFrames/nrm2,AvgColorVec{iFilt},'Color','k','LineWidth',1)
            legend({"Measured EMG" "Measured Avg" FiltLabels{2} FiltLabels{2} FiltLabels{3} FiltLabels{3} "Norm. Force"},'Location','NorthWest')

            subplot(2,1,2)
            nrm1=1000000;
            plot(Time,PW/nrm1)
            plot([Time(1) Time(end)],[TwitchPW TwitchPW]/nrm1)
            legend({FiltLabels{1} FiltLabels{1} FiltLabels{2} FiltLabels{2} FiltLabels{3} FiltLabels{3} "Norm. PW" "Twitch PW"},'Location','NorthWest')

        end                 
    end
end

%% Mean MAV of RC Trials 
clear MAV_Mean
exp_lbl='RC';
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
MeanTime=[8 10]; % Calculating the means at time [8 10]

for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    
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

MAV_Coefs=MAV_Mean./MAV_Mean(1,:);

% MAV_Mean=array2table(MAV_Mean,'VariableNames',[sprintf('PW_%d',),TestFolders(2),TestFolders(3)]);
% Normalized MAV values indicated that the normalizing the EMG signals
% among participants are hard to model. 
% But the sample is 3; the third one might an outlier.

%% Target profile normalization



%% Plot Dropped vs Filtered
clear IndTrials vMVC sMVC IndS IndV DroppedFramesIndexFixed DroppedFrames UnfiltFeat FiltMAV vEMGFeat DroppedFeat DroppedFramesIndexFixed DroppedFrameInd x
close all 
clc

MarginFromDropped=5;  % frames
PlotRange1=[ 5.4 10]; 
PlotRange2=[ 10 15.2]; 
PlotVoli=4; VoliMVCLevels=[10 20 30 40 ];
PlotStim=4; StimMVCLevels=[ 0 10 20 30];
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

    for iFilt=2:length(S.(AnaStruct).AnaPar.FiltLabels)
        FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};
        vEMGLabel=sprintf('%s_vEMG',FeatLabel);
        DroppedFeatLabel=sprintf('%s_vEMG',FeatLabel);

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
            figure(1)
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
                
                plot(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped,ones(length(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped),1)*UnfiltRelatedFramesMean(iPlot,iDropped),...
                    'LineWidth',3,'Color','k')
                hold on
            end

            plot(PlotRangeInd,vEMGFeat(iPlot,:),'.', 'Color','r','LineWidth',0.05)
            plot(PlotRangeInd,UnfiltFeat(iPlot,:),'.','Color','k','LineWidth',0.05)
            plot(DroppedFrames(DroppedFrameInd),DroppedFramesFeat(iPlot,:),'o','Color','b','LineWidth',1.5)
            TrialString=num2str(IndTrials');
            ttl = sprintf('%s filtered vEMG  %s vs Unfiltered at %s, Stim: %d%%, Voli: %d%%, Trials: %s',FiltLabel, FeatLabel,ExpLabel,sMVC,vMVC,TrialString);
            title(ttl);
            legend('Dropped Frame Neighbourhood','Location','northwest')
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'k', 'DisplayName', sprintf('Unfiltered MAV'))
            plot([NaN NaN], [NaN NaN],'.', 'Color', 'r', 'DisplayName', sprintf('%s Filtered MAV',FiltLabel))
            plot([NaN NaN], [NaN NaN],'o', 'Color', 'b', 'DisplayName', sprintf('Dropped Frames'))
            ylabel('EMG MAV');
            xlabel('Frames');
            ylim([0 5*10^-4])
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

        for iPlot=1:length(IndTrials)
            TrialLabel=sprintf('Trial_%d',IndTrials(iPlot));
            
            MeanDrop=mean( PercentDrops);
            StdDrop=std( PercentDrops);

            figure(10)
            subplot(2,1,iFilt-1)
            for iDropped=1:length(DroppedFrames(DroppedFrameInd))
                DroppedInd=DroppedFrames(iDropped);

                plot(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped,ones(length(-MarginFromDropped+DroppedInd:DroppedInd+MarginFromDropped),1)*UnfiltRelatedFramesMean(iPlot,iDropped),...
                'LineWidth',2,'Color','k')
                hold on
                plot(DroppedFrames(DroppedFrameInd),DroppedFramesFeat(iPlot,:),'o','Color','b','LineWidth',1.5)

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
%% Histograms 

close all
for iPlot=1:length(IndTrials)
    for iDropped=1:length(DroppedFrames(DroppedFrameInd))
        DroppedLabel=sprintf('Dropped_%d',iDropped);

        Unfilts(iDropped,:)=UnfiltRelatedFrames.(DroppedLabel){iPlot};
        Filts(iDropped,:)=FiltRelatedFrames.(DroppedLabel){iPlot};
%         Dropped(iPlot)=DroppedFramesFeat(iPlot,:);
    end
    figure
    subplot(3,1,1)
    histogram(Unfilts,10)
    subplot(3,1,2)
    histogram(Filts,10)
    subplot(3,1,3)
    histogram(DroppedFramesFeat,5)
end
  
%% Mean MAV of Occ Trials 

exp_lbl='Occ';
ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
VoliLevels=[10 20 30 40];
StimLevels=[ 0 10 20 30];
MeanTime=[10 15];
% Calculating the means at time [10 15]
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
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

lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
    'PassbandFrequency',5,'PassbandRipple',0.001, ...
    'SampleRate',1/S.Ts);

TrialNum=S.(ExpLabel).iTrial-1;
for iTrial=1:TrialNum
    TrialLabel=sprintf('Trial_%d',iTrial);

    Force=S.(ExpLabel).(TrialLabel).data(:,iForce);

    y1 = filtfilt(lpFilt,Force);
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