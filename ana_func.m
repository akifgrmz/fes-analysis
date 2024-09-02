function S=ana_func(TestFolders,DroppedFiltLabels,NoStimFiltLabels,MAV_MAXMethods,TauTestsForce,TauTestsHand,AvgOcclusionTests,LogModel, BlankTime,gs_orders)
%% mainana_func v1.0


%% Data Inject 

% clc
% clear all
% % TestFolders=["mar20_24"];
% % TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25" "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% % TestFolders=[ "apr20" "oct11" "oct18"];
% % TestFolders=["oct25"];
% % TestFolders=[ "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% % TestFolders=["jan7" "jan11" "jan12"  ];
% % TestFolders=["jun20_24" "jul9_24" "jul21_24"  ]
% % TestFolders=[ "jan7" "jan11" "jan12" "jun20_24" "jul9_24" "jul21_24" "jul31_24" ];
% TestFolders=["feb28_24" "feb29_24" "mar18_24"  "mar20_24" ];
% TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24"];
% TestFolders=[ "jan7" "jan11" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" "jan12" "aug22_24" "aug26_24" "aug29_24"];
% % TestFolders=["aug29_24"]

% for iTest=1:length(TestFolders)
%     TestFiles(iTest)=sprintf("%s_test",TestFolders{iTest});
% end

textt="%s_test";
TestFiles=compose(textt,[TestFolders']);
S = load_test(TestFolders,TestFiles);


Tau_table= readtable('tau_estimates3.csv');
Tau_stats_table= readtable('tau_est_stats.csv');
%% Defining Initial Parameters

NumofTests=length(TestFolders);

BPOrder=ones(1,NumofTests)*20;
fcutlow = ones(1,NumofTests)*10;
fcuthigh = ones(1,NumofTests)*500;

LPPass = ones(1,NumofTests)*20;
LPStop = ones(1,NumofTests)*40;
Ap = ones(1,NumofTests)*.1;

HandLPPass = ones(1,NumofTests)*1;
HandLPStop = ones(1,NumofTests)*100;
HandAp = ones(1,NumofTests)*.1;

lpFiltFpass = ones(1,NumofTests)*2;
lpFiltFstop = ones(1,NumofTests)*10;
lpFiltAp = ones(1,NumofTests)*0.1;
lpFiltAst = ones(1,NumofTests)*30;

FatFiltFpass = ones(1,NumofTests)*0.01;
FatFiltFstop = ones(1,NumofTests)*2;
FatFiltAp = ones(1,NumofTests)*0.1;
FatFiltAst = ones(1,NumofTests)*30;

%% 
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    fs=S.(TestLabel).ExpPar.fs;
    
    % Experiments indices 
    DataIndTable=S.(TestLabel).ExpPar.DataIndTable;
    TableInd=DataIndTable.Properties.VariableNames;
    iForce=table2array(DataIndTable(:,"Force"));
    iTrigger=table2array(DataIndTable(:,"Trigger"));
    iEMG=table2array(DataIndTable(:,"EMG"));
    iTime=table2array(DataIndTable(:,"Time"));
    iPW=table2array(DataIndTable(:,"PW"));

    % Analysis Parameters
    S.(AnaLabel).AnaPar.AnaLabels=["BPFilt_EMG"];
    S.(AnaLabel).AnaPar.AnaInd=table([],'VariableNames',S.(AnaLabel).AnaPar.AnaLabels);
    
    TestNames=strings(1,length(TestFolders));
    if length(char(TestFolders(iTest)))==5 | length(char(TestFolders(iTest)))==4
        TestName=sprintf("%s/23",TestFolders(iTest));
    elseif length(char(TestFolders(iTest)))==8 | length(char(TestFolders(iTest)))==7
        Temp=char(TestFolders(iTest));
        TestName=sprintf("%s/24",Temp(1:5));
    else
        TestName=TestFolders(iTest);
    end

    TestNames(iTest)=TestName;
    S.(AnaLabel).AnaPar.TestName=TestName;

    DataLabels= ["EMG", "Trigger Signal", "Force (N)","Pulse Width","Time (s)"...
        ,"BP Filtered EMG","Blanked EMG","Comb Filtered vEMG","Comb Filtered m-Waves"...
        ,"GS Filtered vEMG","GS Filtered m-Waves"];
    
    FeatLabels=["MAV","MedFreq","MeanFreq","SSC","ZC"];
    FiltLabels=["Unfilt","Comb","GS"];
    S.(AnaLabel).AnaPar.DataLabels=DataLabels;
    S.(AnaLabel).AnaPar.FeatLabels=FeatLabels;
    S.(AnaLabel).AnaPar.FiltLabels=FiltLabels;


    %Incorporate the redo trials at the begginning 
    %Design 10-500Hz 20th order butterworth for filtfilt

    
    % BPOrder=20;
    % fcutlow = 10;
    % fcuthigh = 500;
    
    d1=designfilt('bandpassiir','FilterOrder',BPOrder(iTest), ...
             'HalfPowerFrequency1',fcutlow(iTest),'HalfPowerFrequency2',fcuthigh(iTest), ...
             'SampleRate',fs,'DesignMethod',"butter");
    % fvtool(d1)

    % LPPass = 20;
    % LPStop = 40;
    % Ap = .1;

    d2 = designfilt('lowpassfir','PassbandFrequency',LPPass(iTest),...
      'StopbandFrequency',LPStop(iTest),'PassbandRipple',Ap(iTest),...
      'DesignMethod', 'kaiserwin','SampleRate',fs);
    LPOrder = filtord(d2);

    % HandLPPass = 1;
    % HandLPStop = 100;
    % HandAp = .1;

    d3 = designfilt('lowpassfir','PassbandFrequency',HandLPPass(iTest),...
      'StopbandFrequency',HandLPStop(iTest),'PassbandRipple',HandAp(iTest),...
      'DesignMethod', 'kaiserwin','SampleRate',fs);
    HandLPOrder = filtord(d3);
    % fvtool(d3)
    

    % lpFiltFpass = 2;
    % lpFiltFstop = 10;
    % lpFiltAp = 0.1;
    % lpFiltAst = 30;
    
    stim_freq=S.(TestLabel).ExpPar.FreqList(1);
    lpFilt = designfilt('lowpassfir','PassbandFrequency',lpFiltFpass(iTest),...
      'StopbandFrequency',lpFiltFstop(iTest),'PassbandRipple',lpFiltAp(iTest),...
      'DesignMethod', 'kaiserwin','SampleRate',stim_freq);
    N_lpFilt = filtord(lpFilt);

    % FatFiltFpass = 0.01;
    % FatFiltFstop = 2;
    % FatFiltAp = 0.1;
    % FatFiltAst = 30;
    % 
    FatFilt = designfilt('lowpassfir','PassbandFrequency',FatFiltFpass(iTest),...
      'StopbandFrequency',FatFiltFstop(iTest),'PassbandRipple',FatFiltAp(iTest),...
      'DesignMethod', 'kaiserwin','SampleRate',stim_freq);
    N_FatFilt = filtord(FatFilt);


    
    S.(AnaLabel).AnaPar.BPFilter.fcutlow=fcutlow(iTest);
    S.(AnaLabel).AnaPar.BPFilter.fcuthigh=fcuthigh(iTest);
    S.(AnaLabel).AnaPar.BPFilter.d1=d1;

    S.(AnaLabel).AnaPar.LPFilter.LPPass=LPPass(iTest);
    S.(AnaLabel).AnaPar.LPFilter.LPStop=LPStop(iTest);
    S.(AnaLabel).AnaPar.LPFilter.d2=d2;
    
    S.(AnaLabel).AnaPar.HandFilter.LPPass=HandLPPass(iTest);
    S.(AnaLabel).AnaPar.HandFilter.LPStop=HandLPStop(iTest);

    S.(AnaLabel).AnaPar.lpFilt.lpFiltFpass=lpFiltFpass(iTest);
    S.(AnaLabel).AnaPar.lpFilt.lpFiltFstop=lpFiltFstop(iTest);
    S.(AnaLabel).AnaPar.lpFilt.lpFiltAp=lpFiltAp(iTest);
    S.(AnaLabel).AnaPar.lpFilt.lpFiltAst=lpFiltAst(iTest);

    S.(AnaLabel).AnaPar.FatFilt.FatFiltFpass=FatFiltFpass(iTest);
    S.(AnaLabel).AnaPar.FatFilt.FatFiltFstop=FatFiltFstop(iTest);
    S.(AnaLabel).AnaPar.FatFilt.FatFiltAp=FatFiltAp(iTest);
    S.(AnaLabel).AnaPar.FatFilt.FatFiltAst=FatFiltAst(iTest);


    S.(AnaLabel).AnaPar.FatFilt.FatFilt=FatFilt;
    S.(AnaLabel).AnaPar.FatFilt.N_FatFilt=N_FatFilt;
    S.(AnaLabel).AnaPar.lpFilt.lpFilt=lpFilt;
    S.(AnaLabel).AnaPar.LPFilter.LPOrder=LPOrder;
    S.(AnaLabel).AnaPar.HandFilter.LPOrder=HandLPOrder;
    S.(AnaLabel).AnaPar.BPFilter.BPOrder=BPOrder;    
    S.(AnaLabel).AnaPar.lpFilt.N_lpFilt=N_lpFilt;

    S.(AnaLabel).AnaPar.LPFilters{1}=d2;
    S.(AnaLabel).AnaPar.LPFilters{2}=d3;

    % Target Line Generation
    Target=[ linspace(0,0,ceil(5*fs)+1) linspace(0,1,ceil(5*fs)) linspace(1,1,ceil(7*fs)-2)]';
    length(Target);
    ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).('Occ'));
    S.(AnaLabel).(ExpLabel).Target=Target;
end

%
%% Fixing PW invalid values issue (-inf)
%

TrigThres=1;
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    S.(TestLabel).ExpPar.TrigThres=TrigThres;
    iPW=S.(TestLabel).ExpPar.DataIndTable.("PW");
    iTrigger=S.(TestLabel).ExpPar.DataIndTable.("Trigger");
    ExpRuns=S.(TestLabel).ExpRuns;
    EffortType=S.(TestLabel).ExpPar.EffortType;

    for iExp=1:length(S.(TestLabel).ExpPar.ExpLabels)
        ExpLabel=S.(TestLabel).ExpPar.ExpLabels(iExp);
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;

        if ExpRuns(iExp)
            RedoTrials=S.(TestLabel).(ExpLabel).RedoTrials;
            if ~isempty(RedoTrials)
                for iRedo=1:length(RedoTrials)
                    TrialLabel=sprintf('Trial_%d',RedoTrials(iRedo));
                    RedoLabel=sprintf('RedoTrial_%d',RedoTrials(iRedo));

                    S.(TestLabel).(ExpLabel).(TrialLabel).data=S.(TestLabel).(ExpLabel).(RedoLabel).data; %% Incorporate redo trials
                    S.(AnaLabel).(ExpLabel).(TrialLabel).Redo=true;
                    
                    x=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,iPW);
                    y=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,iTrigger);
                        % # PW
                    x(x<0)=0;
                    x(isnan(x))=0;
                    StimOnInd(x~=0)=true;

                    S.(TestLabel).(ExpLabel).(TrialLabel).data(:,iPW)=x;
                    S.(AnaLabel).(ExpLabel).(TrialLabel).StimInd=StimOnInd;
                    StimOnInd=[];
                                % # trigger
                    y(y>TrigThres)=1;
                    y(y<TrigThres)=0;
                    S.(TestLabel).(ExpLabel).(TrialLabel).data(:,iTrigger)=y;

%                     dt=S.(TestStruct).(ExpLabel).(TrialLabel).data;
%                     [r,~]=size(dt);
%                     T=array2table(dt,'VariableNames',S.(TestStruct).ExpPar.DataInd);
% 
%                     S.(TestStruct).(ExpLabel).(TrialLabel).data=T;
                end
            end 
            
            if ExpLabel=="OccTrials"
                NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
            end
        
            for iTrial=1:NumofTrials
                TrialLabel=sprintf('Trial_%d',iTrial);

                S.(AnaLabel).(ExpLabel).(TrialLabel).data=S.(TestLabel).(ExpLabel).(TrialLabel).data;
                
                x=S.(AnaLabel).(ExpLabel).(TrialLabel).data(:,iPW);
                y=S.(AnaLabel).(ExpLabel).(TrialLabel).data(:,iTrigger);
                    % # PW
                x(x<0)=0;
                x(isnan(x))=0;
                StimOnInd(x~=0)=true;

                S.(AnaLabel).(ExpLabel).(TrialLabel).data(:,iPW)=x;
                S.(AnaLabel).(ExpLabel).(TrialLabel).StimInd=StimOnInd;
                StimOnInd=[];

                            % # trigger
                y(y>TrigThres)=1;
                y(y<TrigThres)=0;
                S.(AnaLabel).(ExpLabel).(TrialLabel).data(:,iTrigger)=y;
        
                %% Making "data" a table
                    data=S.(AnaLabel).(ExpLabel).(TrialLabel).data;
                    [r,c]=size(data);
                    T=array2table(data,'VariableNames',S.(TestLabel).ExpPar.DataInd(:,1:c));

                    if EffortType=="Force"
                        T.Hand=zeros(r,1)*NaN;
                    end
                 
                S.(AnaLabel).(ExpLabel).(TrialLabel).data=T;
                if ~isfield(S.(AnaLabel).(ExpLabel).(TrialLabel), 'Redo')
                    S.(AnaLabel).(ExpLabel).(TrialLabel).Redo=false;
                end
            end
        end
    end
end

%%EMG and Force Filtering 
% ExpLabels=S.(TestStruct).ExpPar.ExpLabels.()

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    AnaLabel=sprintf("%s_ana",TestFolders(iTest));
    AmpGain=S.(TestLabel).ExpPar.AmpGain;
    EffortLabel=S.(TestLabel).ExpPar.EffortType;

    ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    ExpRuns=S.(TestLabel).ExpRuns;
    ExpRuns(3)=false;  %% No Custom Trials
    ExpstoAna=ExpLabels(ExpRuns);
    
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna(iExp);

        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        % NumofTrials(ExpLabel=="OccTrials")=S.(TestStruct).(ExpLabel).ListedNumofTrials;
        
        if ExpLabel=="OccTrials"
            NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
  
            EMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('EMG')/AmpGain;
            EffortMeasure=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortLabel);
            LP=S.(AnaLabel).AnaPar.LPFilter.d2;
            BP=S.(AnaLabel).AnaPar.BPFilter.d1;

            BPFilt_EMG = filtfilt(BP,EMG); % BP filtering for EMG
            Filt_EffortMeasure = filtfilt(LP,EffortMeasure);% LP Filtering for force 

            S.(AnaLabel).(ExpLabel).(TrialLabel).data.('BPFilt_EMG')=BPFilt_EMG;
            S.(AnaLabel).(ExpLabel).(TrialLabel).data.(sprintf('Filt_%s',EffortLabel))=Filt_EffortMeasure;
        end
    end
end

%%Trigger and Blanking

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    ExpRuns=S.(TestLabel).ExpRuns;
    ExpRuns(3)=false;  %% No Custom Trials
    ExpstoAna=ExpLabels(ExpRuns);
    
    BlankLength=round(BlankTime(iTest)*fs)+1;
    S.(AnaLabel).AnaPar.BlankLength=BlankLength;
    S.(AnaLabel).AnaPar.BlankTime=BlankTime;
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna(iExp);
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;

        if ExpLabel=="OccTrials"
            NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);

            y=S.(AnaLabel).(ExpLabel).(TrialLabel).data.("BPFilt_EMG");
            x=S.(AnaLabel).(ExpLabel).(TrialLabel).data.("Trigger");

            Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.("Time"); 
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
    
            S.(AnaLabel).(ExpLabel).(TrialLabel).Indices=table([FallingInd], [RisingInd]);
            S.(AnaLabel).(ExpLabel).(TrialLabel).data.("BlankedEMG")=z;
            S.(AnaLabel).(ExpLabel).(TrialLabel).data.("Time")=Time;
            S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames=BegofFrames;
            S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength=FrameLength;

        end
    end
end

%% Plots after filtering and preproccess 
% clc
% close all
% Trials=[ 11 12 13];  % Trial
% TimeRange=[1 12];  % in seconds
% % TestLabel=sprintf('%s_test',FolderNames);
% ylims=5;
% ExpNum=4;
% for iTest=1:length(TestFolders)
%     TestLabel=sprintf('%s_test',TestFolders(iTest));
%     AnaLabel=sprintf('%s_ana',TestFolders(iTest));
% 
%     TestName=TestFolders(iTest);
%     ExpLabel=S.(TestLabel).ExpPar.ExpLabels(ExpNum);
%     EffortType=S.(TestLabel).ExpPar.EffortType;
% 
%     for iTrial=1:length(Trials)
%         TrialLabel=sprintf('Trial_%d',Trials(iTrial));
%         DataIndTable= S.(TestLabel).ExpPar.DataIndTable;
% 
%         TimeRange=S.(TestLabel).(ExpLabel).StimRange;
%         Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
%         TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
% 
%         % RawEMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('EMG');
%         BlankedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('BlankedEMG');
% 
%         Trigger=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('Trigger');
%         EffortMeas=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).(EffortType);
%         PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('PW');
% 
%         figure(iTest)
%         subplot(length(Trials),1,iTrial)
%         plot(Time(TimeInd),Trigger/2000,'r','LineWidth',1)
%         hold
%         plot(Time(TimeInd),BlankedEMG,'b','LineWidth',2)
%         legend({'Trigger(a.u.)','Raw EMG (mV)'})
%         title(sprintf('Test:%s, Trial: %d', TestName, Trials(iTrial)))
%         xlabel('Time (s)')
%         % ylim([-1 1]*10^-3)
%     end
% OccTable=S.(TestLabel).OccTrials.RepTableMat;
% OccTable=[[1:length(OccTable)]' OccTable]
% 
% end

%% Frames Matrix for Force and EMG
% (also remove blanked periods)
% (Add unfilt as a filter and add e m-wave as zero matrix)
% Identify Zero stim trials -> Done

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    EffortType=S.(TestLabel).ExpPar.EffortType;

    ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    ExpRuns=S.(TestLabel).ExpRuns;
    ExpRuns(3)=false;  %% No Custom Trials
    ExpstoAna=ExpLabels(ExpRuns);

    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            
            PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('PW');
            BegofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);

            clear PWofFrames
            for iFrame=1:NumofFrames
                %PW
                PWofFrames(iFrame)=PW(BegofFrames(iFrame));
            end

            S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames=PWofFrames';
        end
    end

%     ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
%     ExpstoAna=ExpLabels([1,2,3,4,5]);
    BlankLength=S.(AnaLabel).AnaPar.BlankLength;
    
    for iExp=1:length(ExpstoAna)
        ExpLabel=ExpstoAna{iExp};
        
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;

        if ExpLabel=="OccTrials"
            NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
        
        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            FrameLength=S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength;
            
            BlankEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('BlankedEMG');
            BPFilt_EMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('BPFilt_EMG');
            
            Force=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Force');
            Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
            Trigger=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Trigger');

            Hand=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Hand');

            BegofFrames= S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);

            clear y x f t tg FrameLengths PWFrames TmFrames h
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
                % # Trigger Frames
                tg(:,iFrame)= Trigger(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
                    % # Hand
                h(:,iFrame)= Hand(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
            end

            S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames=y;
            S.(AnaLabel).(ExpLabel).(TrialLabel).EMGFrames=x;
            S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames=f;
            S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames=t;
            S.(AnaLabel).(ExpLabel).(TrialLabel).TriggerFrames=tg;
            S.(AnaLabel).(ExpLabel).(TrialLabel).BlankTmFrames=t;
            S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLengths=FrameLengths;

            S.(AnaLabel).(ExpLabel).(TrialLabel).HandFrames=h;
        end
    end
    
    % Target Frames
    ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).('Occ'));
    Target=S.(AnaLabel).(ExpLabel).Target;    
    BegofFrames= S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
    NumofFrames=length(BegofFrames);
    TargetFrames=[];
    for iFrame=1:NumofFrames-1
        % # Target Line Frames
        TargetFrames(:,iFrame)=Target(BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1));
    end
    S.(AnaLabel).(ExpLabel).TargetFrames=TargetFrames;
end


%%Extracting Dropped Frames 

ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ'); 
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});

    StimRange=S.(TestLabel).(ExpLabel).StimRange;

    if ~S.(TestLabel).(ExpLabel).dropped 
        continue;
    else

        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        stim_freq=S.(TestLabel).ExpPar.FreqList(1);
        TrialsPW=S.(TestLabel).(ExpLabel).TrialsPW;
        TrialsPW=S.(TestLabel).(ExpLabel).RepTableMat(:,5); % 5--> PW

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);

            % StimRange=S.(TestLabel).(ExpLabel).(TrialLabel).StimRange;

            FrameRange=[(stim_freq+1)*StimRange(1),stim_freq*StimRange(2)];
            PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('PW');

            BegofFrames= S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
            NumofFrames=length(BegofFrames);

            if TrialsPW(iTrial) == 0
                ZeroInd=zeros(NumofFrames,1);
            else
                ZeroInd=find(( PW(BegofFrames)==0)==1);
            end

            DroppedFrames= ZeroInd(ZeroInd>=FrameRange(1) & ZeroInd<=FrameRange(2) );
            DroppedEMG= S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames(:,DroppedFrames);

            S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames=DroppedFrames;
            S.(AnaLabel).(ExpLabel).(TrialLabel).ZeroInd=ZeroInd;
            S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG=DroppedEMG;
        end
    end
end

%% Plotting dropped frames
% 
% S = load_test ; 
% %%
% close all
% Lbl='Occ';
% PlotTrial=[ 15 15];
% PlotFrame=floor([ 488 488]);
% 
% for iTest=1:length(TestFolders)
%     TestStruct=sprintf("%s_test",TestFolders{iTest});
%     AnaStruct=sprintf("%s_ana",TestFolders{iTest});
% 
%     ExpTable=S.(AnaStruct).AnaPar.ExpTable.(Lbl);
%     ExpLabel=ExpTable{1};
%     DataLabels=S.(AnaStruct).AnaPar.DataLabels;
%     AmpGain=S.(TestStruct).ExpPar.AmpGain;
% 
%     for iTrial=PlotTrial(1):PlotTrial(2)
%         TrialLabel=sprintf('Trial_%d',iTrial);
%         DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
% 
%         if isempty(DroppedFrames)
%             subplot(length(TestFolders),1,iTest)
%             ttl=sprintf('BP Filtered %s at %s,Test: %s, TrialNum: %d, Frame: %d-%d,'...
%             ,DataLabels{1},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2));
%             title(ttl);
%             text(0.3,0.5,'No Dropped Frames for This Trial')
% 
%             continue 
%         end
% 
%         FrameLength=S.(AnaStruct).(ExpLabel).(TrialLabel).FrameLength;
%         BegofFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).BegofFrames;
%     %     Ind=[BlankLength+PlotFrame(1)*FrameLength BlankLength+(PlotFrame(2))*FrameLength];
%         Ind= [BegofFrames(PlotFrame(1)) BegofFrames(PlotFrame(2)+1)];
% 
%         BPFilt_EMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BPFilt_EMG');
%         BlankEMG=S.(AnaStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('BlankedEMG');
%         Trig=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Trigger');
%         EMG=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('EMG');
%         Time=S.(TestStruct).(ExpLabel).(TrialLabel).data(Ind(1):Ind(2),:).('Time');
% 
%         figure(1)
%         subplot(length(TestFolders),1,iTest)
%         plot(Time,Trig/1000,'b','LineWidth',2)
%         hold on
%         plot(Time,EMG/AmpGain,'k','LineWidth',2)
%         plot(Time,BlankEMG,'r','LineWidth',2)
%         legend({'Trigger','Before Blanking and Filtering', 'After Blanking and Filtering'})
%         ttl=sprintf('BP Filtered %s at %s,Test: %s, TrialNum: %d, Frame: %d-%d,'...
%             ,DataLabels{1},ExpLabel,TestFolders{iTest},iTrial,PlotFrame(1),PlotFrame(2));
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
%         ylim([1.5*min(BPFilt_EMG) 1.5*max(BPFilt_EMG)])
%     end
% end



%% M-wave filtering  
% (Remove blanked periods before filtering, Do unfiltered Frames)
% (update to recreate the vector forms of the signals)

ExpstoAna= ["MVC","RC","Occ","Fat","Ramp"];
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    
    FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
    gs_order=gs_orders(iTest);

    for iExp=1:length(ExpstoAna)
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(ExpstoAna(iExp));
        
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;

        %%%%%%%%% ------------------------------------ must go to tidy data 
        if ExpLabel=="OccTrials"
            NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
        
        if isfield(S.(TestLabel).(ExpLabel), 'dropped')
            
            dropped=S.(TestLabel).(ExpLabel).dropped;
        else
            dropped=0;
        end

        if isfield(S.(TestLabel).(ExpLabel), 'TrialsPW')
            
            TrialsPW=S.(TestLabel).(ExpLabel).TrialsPW;
        else
            TrialsPW(1:NumofTrials)=0; 
        end

        if  ExpLabel==S.(TestLabel).ExpPar.ExpTable(1,:).('Occ')

            TrialsPW=S.(TestLabel).(ExpLabel).RepTableMat(:,5);
        end

         %%%%%%%%% -------------------------------------- must go to tidy data 

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            clear T Time TimeFrames x_frames
            BegofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames;
            FrameNum=length(BegofFrames);
            
            if dropped==1 && TrialsPW(iTrial) ~= 0
                DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
                KeepInd=setdiff(1:FrameNum-1,DroppedFrames);
                KeepInd=KeepInd(1:end-2);
                x_frames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames(:,KeepInd); % excluding dropped frames 
                [FrameLength, FrameNum]=size(x_frames);
                TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankTmFrames(:,KeepInd);
                Time=reshape(TimeFrames,[1],FrameLength*FrameNum);
                S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength=FrameLength;
            
            elseif dropped==1 && TrialsPW(iTrial) == 0
                DroppedFrames=S.(AnaLabel).(ExpLabel).Trial_1.DroppedFrames;
                KeepInd=setdiff(1:FrameNum-1,DroppedFrames);
                
                x_frames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames(:,KeepInd); % excluding dropped frames 
                [FrameLength, FrameNum]=size(x_frames);
                TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankTmFrames(:,KeepInd);
                Time=reshape(TimeFrames,[1],FrameLength*FrameNum);
                S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength=FrameLength;
                
            else 
                x_frames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames;
                [FrameLength, FrameNum]=size(x_frames);

                KeepInd=[];
                TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames;
                Time=reshape(TimeFrames,[1],FrameLength*FrameNum);
            end

            x_frames_withdropped=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames; 
            [FrameLength, FrameNum_withdropped]=size(x_frames_withdropped);

            S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd=KeepInd;
            S.(AnaLabel).(ExpLabel).Time=table(Time','VariableName',"Time");

            % # Unfilt
            iFilt=1;
            FiltLabel=FiltLabels(iFilt);
            UnfiltvEMGFrames=x_frames;
            UnfiltMWaveFrames=zeros(size(UnfiltvEMGFrames));
            UnfiltvEMG=reshape(UnfiltvEMGFrames,1,FrameLength*FrameNum);
            UnfiltMWave=reshape(UnfiltMWaveFrames,1,FrameLength*FrameNum);

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=UnfiltMWaveFrames;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=UnfiltvEMGFrames;

            % With dropped filtered
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped=zeros(size(x_frames_withdropped));
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped=x_frames_withdropped;            
            
%             T=[table(UnfiltvEMG','VariableName',"UnfiltvEMG") table(UnfiltMWave','VariableName',"UnfiltMWave")];
%             S.(AnaLabel).(ExpLabel).(FiltLabel).T=T;   


            if dropped==1 && TrialsPW(iTrial) ~= 0
                BlankEMGFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=UnfiltvEMGFrames;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=UnfiltMWaveFrames;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped=BlankEMGFrames;
            end

            % # Comb Filter
            iFilt=2;

            [CombMWaveFrames,CombvEMGFrames]=FiltComb(x_frames);
            CombvEMG=reshape(CombvEMGFrames,1,FrameLength*FrameNum);
            CombMWave=reshape(CombMWaveFrames,1,FrameLength*FrameNum);

            FiltLabel=FiltLabels{iFilt};
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=CombMWaveFrames;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=CombvEMGFrames;
%             T=[table(CombvEMG','VariableName',"CombvEMG") table(CombMWave','VariableName',"CombMWave")];
%             S.(AnaLabel).(ExpLabel).(FiltLabel).T=T;   

            % Dropped frames filtered
            [CombMWaveFrames_filtdropped,CombvEMGFrames_filtdropped]=FiltComb(x_frames_withdropped);

            CombvEMG_filtdropped=reshape(CombvEMGFrames_filtdropped,1,FrameLength*FrameNum_withdropped);
            CombMWave_filtdropped=reshape(CombMWaveFrames_filtdropped,1,FrameLength*FrameNum_withdropped);

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped=CombMWaveFrames_filtdropped;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped=CombvEMGFrames_filtdropped;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWave_filtdropped=CombMWave_filtdropped;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMG_filtdropped=CombvEMG_filtdropped;

            if dropped==1 && TrialsPW(iTrial) ~= 0

                BlankEMGFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=CombvEMGFrames;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=CombMWaveFrames;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped=BlankEMGFrames;

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

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames=vEMGMat;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames=MWaveMat;

            GSvEMG=reshape(vEMGMat,1,FrameLength*FrameNum);
            GSMWave=reshape(MWaveMat,1,FrameLength*FrameNum);

%             T=[table(GSvEMG','VariableName',"GSvEMG") table(GSMWave','VariableName',"GSMWave")];
%             S.(AnaLabel).(ExpLabel).(FiltLabel).T=T;   

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMG=GSvEMG;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWave=GSMWave;

            if dropped==1 && TrialsPW(iTrial) ~= 0

                BlankEMGFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=vEMGMat;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped=BlankEMGFrames;
                BlankEMGFrames(:,KeepInd)=MWaveMat;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped=BlankEMGFrames;
            end

            % Dropped frames filtered
            clear vEMGMat_withdropped MWaveMat_withdropped
            for iFrame = 1:FrameNum_withdropped-gs_order
                temp_vect = SUB_GS_filter(x_frames_withdropped(:,iFrame:gs_order+iFrame),FrameLength,gs_order);
                vEMGMat_withdropped(:,iFrame+gs_order) = temp_vect(:,1);
                MWaveMat_withdropped(:,iFrame+gs_order) = temp_vect(:,2);
            end

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped=vEMGMat_withdropped;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped=MWaveMat_withdropped;

            GSvEMG_filtdropped=reshape(vEMGMat_withdropped,1,FrameLength*FrameNum_withdropped);
            GSMWave_filtdropped=reshape(MWaveMat_withdropped,1,FrameLength*FrameNum_withdropped);

            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMG_filtdropped=GSvEMG_filtdropped;
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWave_filtdropped=GSMWave_filtdropped;


            % # Blanking filter comes here
        end
    end
    S.(AnaLabel).AnaPar.gs_order=gs_order;
end


% %% Plotting, Debugging 
% % close all
% Lbl='Occ';
% PlotTrial=[ 10 10];
% PlotTime=[ 12 13];
% % PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);
% PlotFrame=floor([ 478 482]);
% 
% TeststoPlot=TestFolders();
% 
% for iTest=1:length(TeststoPlot)
%     TestLabel=sprintf("%s_test",TeststoPlot(iTest));
%     AnaLabel=sprintf("%s_ana",TeststoPlot(iTest));
%     DataLabels=S.(AnaLabel).AnaPar.DataLabels;
%     TestName=S.(AnaLabel).AnaPar.TestName;
%     ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).(Lbl));
% 
%     % PlotFrame=floor([ 300 305]); % first frame is skipped
%     % ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
%     % ExpstoAna=ExpLabels([2,3,4,5]);
%     FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
% 
%     clear EMGFrames MWaveFrames vEMG MWave MWaveFrames_filtdropped vEMGFrames_filtdropped vEMG_filtdropped MWave_filtdropped
%     % vEMGFrames_filtdropped=zeros(FrameLength,PlotFrame(2)-PlotFrame(1)+1,length(FiltLabels));
%     for iTrial=PlotTrial(1):PlotTrial(2)
%         TrialLabel=sprintf('Trial_%d',iTrial);
%         FrameLength=S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength;
% 
%         for iFilt=1:length(FiltLabels)
%             FiltLabel=FiltLabels{iFilt};
% 
%             EMGFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
%             MWaveFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
% 
%             vEMG(:,iFilt)=reshape(EMGFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
%             MWave(:,iFilt)=reshape(MWaveFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
% 
%             vEMGFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
%             MWaveFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
% 
%             vEMG_filtdropped(:,iFilt)=reshape(vEMGFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
%             MWave_filtdropped(:,iFilt)=reshape(MWaveFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
% 
%             % vEMGFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
%             % MWaveFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
%             % 
%             % vEMG_withdropped(:,iFilt)=reshape(vEMGFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
%             % MWave_withdropped(:,iFilt)=reshape(MWaveFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
%         end
% 
%         TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
%         TriggerFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TriggerFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
% 
%         Time=reshape(TimeFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
%         Trig=reshape(TriggerFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
% 
%         figure('Name',sprintf('%s Results, Test: %s',FiltLabels{2},TeststoPlot(iTest)),'NumberTitle','off')
%         subplot(2,1,1)
%         plot(Time,Trig/10000,'b','LineWidth',2)
%         hold on
%         plot(Time,vEMG(:,1),'k','LineWidth',2)
% 
%         plot(Time,vEMG(:,2),'r','LineWidth',2)
%         legend({'Trigger(a.u)', 'Unfiltered', 'Comb vEMG'})
%         ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d'...
%             ,DataLabels{8},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%         subplot(2,1,2)
%         plot(Time,vEMG(:,1),'k','LineWidth',2)
%         hold
%         plot(Time,MWave(:,2),'r','LineWidth',2)
%         legend({'Unfiltered', 'Comb MWaves'})
%         ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%             DataLabels{9},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%         figure('Name',sprintf('%s Results, Test: %s',FiltLabels{3},TeststoPlot(iTest)),'NumberTitle','off')
%         subplot(2,1,1)
%         plot(Time,Trig/10000,'b','LineWidth',2)
%         hold on
%         plot(Time,vEMG(:,1),'k','LineWidth',2)
% 
%         plot(Time,vEMG(:,3),'r','LineWidth',2)
%         legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
%         ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%             DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%         subplot(2,1,2)
%         plot(Time,vEMG(:,1),'k','LineWidth',2)
%         hold
%         plot(Time,MWave(:,3),'r','LineWidth',2)
%         legend({'Unfiltered', 'GS MWaves'})
%         ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%             DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%         figure(10)
%         subplot(2,1,1)
%         plot(Time,Trig/10000,'b','LineWidth',2)
%         hold on
%         plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
% 
%         plot(Time,vEMG_filtdropped(:,3),'r','LineWidth',2)
%         legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
%         ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%             DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%         subplot(2,1,2)
%         plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
%         hold
%         plot(Time,MWave_filtdropped(:,3),'r','LineWidth',2)
%         legend({'Unfiltered', 'GS MWaves'})
%         ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
%             DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
%         title(ttl);
%         xlabel(DataLabels{5})
%         ylabel(DataLabels{1})
% 
%     end
% end

% %% Plotting the dropped frames 
% close all
% iTrial=16;
% TrialLabel=sprintf('Trial_%d',iTrial);
% ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
% FiltLabel="GS";
% for iTest=1:length(TestFolders)
%     TestStruct=sprintf("%s_test",TestFolders{iTest});
%     AnaLabel=sprintf("%s_ana",TestFolders{iTest});
%     if ~S.(TestStruct).(ExpLabel).dropped 
%                     figure(iTrial)
%             subplot(length(TestFolders),1,iTest)
%             text(0.3,0.5,'No Dropped Frames for This Trial')
%             ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
%             title(ttl);
%         continue;
%     else
%        DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;   
%         if isempty(DroppedFrames)
%             figure(iTrial)
%             subplot(length(TestFolders),1,iTest)
%             text(0.3,0.5,'No Dropped Frames for This Trial')
%             ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
%             title(ttl);
%             continue;
%         end
% 
%         DataLabels=S.(AnaLabel).AnaPar.DataLabels;
%         DroppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG(1:end-3,:);
%         FiltDropped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(:,DroppedFrames);
%         clear lgd
% 
%         for iFrame=1:length(DroppedFrames)
%             figure(iTrial)
%             subplot(length(TestFolders),2,iTest*2-1)
%             plot(DroppedEMG(:,iFrame),'LineWidth',2)
%             hold on
%             lgd(iFrame)=sprintf("Frame %d",DroppedFrames(iFrame));
%             legend(lgd,'Location','NorthWest');
%             ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
%             title(ttl);
%             xlabel('Samples')
%             ylabel('BP Filtered EMG')
%             subplot(length(TestFolders),2,iTest*2)
%             plot(FiltDropped(:,iFrame),'LineWidth',2)
%             hold on 
%             ttl=sprintf('Filtered Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
%             title(ttl);
%             ylabel('BP Filtered and m-Wave Filtered')
% 
% 
%             % ylim([3 -3]*10^-3)
% 
%         end
%     end
% end


%% EMG Features
% (move force averaging to the top where BP is held)
% lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
%     'PassbandFrequency',1,'PassbandRipple',0.1,'SampleRate',stim_freq);

AnaLabel=sprintf("%s_ana",TestFolders(1));

clear FiltFramesInd ExpLabels
ExpstoAna= ["MVC";"RC" ;"Occ";"Fat";"Ramp"];

FieldNames=["MWaveFrames", "MWaveFrames_filtdropped";
            "vEMGFrames","vEMGFrames_filtdropped"];  %fields that contain EMG that get analized here
FieldSuffix=["","_filtdropped"];

for iExp=1:length(ExpstoAna)
    ExpLabels(iExp)=S.(TestLabel).ExpPar.ExpTable(1,:).(ExpstoAna(iExp));
end

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    lpFilt=S.(AnaLabel).AnaPar.lpFilt.lpFilt;
    FatFilt=S.(AnaLabel).AnaPar.FatFilt.FatFilt;

    fs=S.(TestLabel).ExpPar.fs;
    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels(iExp);
        
        TrialNum=S.(TestLabel).(ExpLabel).NumofTrials;
        if ExpLabel=="OccTrials"
            TrialNum=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
        
        for iTrial=1:TrialNum
            TrialLabel=sprintf('Trial_%d',iTrial);

            x1Frames=S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames;

            FiltForceFrames = filtfilt(lpFilt,mean(x1Frames)); % mean force for each frame so fs=35hz
            S.(AnaLabel).(ExpLabel).(TrialLabel).FiltForceFrames=FiltForceFrames;

            for iFilt=1:length(S.(AnaLabel).AnaPar.FiltLabels)
                FiltLabel=S.(AnaLabel).AnaPar.FiltLabels{iFilt};
                
                for iField=1:length(FieldNames)
                    FieldLabel=FieldNames(:,iField);

                    vEMGFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(FieldLabel(2));
                    MWaveFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(FieldLabel(1));
                    [~,FrameNum]=size(vEMGFrames);
    
                    MAV_vEMG=mean(abs(vEMGFrames)); 
                    MAV_MWaves=mean(abs(MWaveFrames)); 
    
                    MedFreq_vEMG=zeros(1,FrameNum);
                    MeanFreq_vEMG=zeros(1,FrameNum);
                    MedFreq_MWaves=zeros(1,FrameNum);
                    MeanFreq_MWaves=zeros(1,FrameNum);
    
                    SSC_vEMG=zeros(1,FrameNum);
                    ZC_vEMG=zeros(1,FrameNum);
                    SSC_MWave=zeros(1,FrameNum);
                    ZC_MWave=zeros(1,FrameNum);
    
                    for iFrame=1:FrameNum
    
                        [MedFreq_vEMG(iFrame), MeanFreq_vEMG(iFrame)]=MedMeanFreq(vEMGFrames(:,iFrame),fs);
                        [MedFreq_MWaves(iFrame), MeanFreq_MWaves(iFrame)]=MedMeanFreq(MWaveFrames(:,iFrame),fs);
    
                        SSC_vEMG(iFrame)=NumSsc(vEMGFrames(:,iFrame));
                        SSC_MWave(iFrame)=NumSsc(MWaveFrames(:,iFrame));
    
                        ZC_vEMG(iFrame)=NumZc(vEMGFrames(:,iFrame));
                        ZC_MWave(iFrame)=NumZc(MWaveFrames(:,iFrame));
    
                    end
    
                    MAV_vEMG(isnan(MAV_vEMG))=0;
                    MAV_MWaves(isnan(MAV_MWaves))=0;
    
                    MedFreq_vEMG(isnan(MedFreq_vEMG))=0;
                    MeanFreq_vEMG(isnan(MeanFreq_vEMG))=0;
                    MedFreq_MWaves(isnan(MedFreq_MWaves))=0;
                    MeanFreq_MWaves(isnan(MeanFreq_MWaves))=0;
                    SSC_vEMG(isnan(SSC_vEMG))=0;
                    SSC_MWave(isnan(SSC_MWave))=0;
                    ZC_vEMG(isnan(ZC_vEMG))=0;
                    ZC_MWave(isnan(ZC_MWave))=0;
    
                    Label=sprintf("Feats%s",FieldSuffix(iField));
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                         table(MAV_vEMG',MedFreq_vEMG',MeanFreq_vEMG',SSC_vEMG',...
                         ZC_vEMG',MAV_MWaves',MedFreq_MWaves',MeanFreq_MWaves',...
                         SSC_MWave',ZC_MWave','VariableNames',["MAV_vEMG" "MedFreq_vEMG"...
                         "MeanFreq_vEMG" "SSC_vEMG" "ZC_vEMG" "MAV_MWave"...
                         "MedFreq_MWave" "MeanFreq_MWave" "SSC_MWave" "ZC_MWave"]);
                    
                     if iExp ~= 1 % dont do this for MVC for having such small trial times less than the filter length
                         
                        FiltMAV_vEMG = filtfilt(lpFilt,MAV_vEMG);
                        FiltMedFreq_vEMG = filtfilt(lpFilt,MedFreq_vEMG);
                        FiltMeanFreq_vEMG = filtfilt(lpFilt,MeanFreq_vEMG);
                        FiltMAV_MWaves = filtfilt(lpFilt,MAV_MWaves);
                        FiltMedFreq_MWaves = filtfilt(lpFilt,MedFreq_MWaves);
                        FiltMeanFreq_MWaves = filtfilt(lpFilt,MeanFreq_MWaves);
                        FiltSSC_vEMG = filtfilt(lpFilt,SSC_vEMG);
                        FiltZC_vEMG = filtfilt(lpFilt,ZC_vEMG);
                        FiltSSC_MWave = filtfilt(lpFilt,SSC_MWave);
                        FiltZC_MWave = filtfilt(lpFilt,ZC_MWave);
    
                        Label=sprintf("FiltFeats%s",FieldSuffix(iField));
                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                            table(FiltMAV_vEMG',FiltMedFreq_vEMG',FiltMeanFreq_vEMG',...
                            FiltSSC_vEMG',FiltZC_vEMG',FiltMAV_MWaves',FiltMedFreq_MWaves',...
                            FiltMeanFreq_MWaves',FiltSSC_MWave',FiltZC_MWave','VariableNames',...
                            [ "Filt_MAV_vEMG" "Filt_MedFreq_vEMG" "Filt_MeanFreq_vEMG"...
                            "Filt_SSC_vEMG" "Filt_ZC_vEMG" "Filt_MAV_MWave" "Filt_MedFreq_MWave"...
                            "Filt_MeanFreq_MWave" "Filt_SSC_MWave" "Filt_ZC_MWave"]);
    
                        FatMAV_vEMG = filtfilt(FatFilt,MAV_vEMG);
                        FatMedFreq_vEMG = filtfilt(FatFilt,MedFreq_vEMG);
                        FatMeanFreq_vEMG = filtfilt(FatFilt,MeanFreq_vEMG);
                        FatMAV_MWaves = filtfilt(FatFilt,MAV_MWaves);
                        FatMedFreq_MWaves = filtfilt(FatFilt,MedFreq_MWaves);
                        FatMeanFreq_MWaves = filtfilt(FatFilt,MeanFreq_MWaves);
                        FatSSC_vEMG = filtfilt(FatFilt,SSC_vEMG);
                        FatZC_vEMG = filtfilt(FatFilt,ZC_vEMG);
                        FatSSC_MWave = filtfilt(FatFilt,SSC_MWave);
                        FatZC_MWave = filtfilt(FatFilt,ZC_MWave);
    
                        Label=sprintf("FatFeats%s",FieldSuffix(iField));                    
                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                            table(FatMAV_vEMG',FatMedFreq_vEMG',FatMeanFreq_vEMG',FatSSC_vEMG',...
                            FatZC_vEMG',FatMAV_MWaves',FatMedFreq_MWaves',...
                            FatMeanFreq_MWaves',FatSSC_MWave',FatZC_MWave','VariableNames',...
                            [ "Filt_MAV_vEMG" "Filt_MedFreq_vEMG" "Filt_MeanFreq_vEMG"...
                            "Filt_SSC_vEMG" "Filt_ZC_vEMG" "Filt_MAV_MWave" "Filt_MedFreq_MWave"...
                            "Filt_MeanFreq_MWave" "Filt_SSC_MWave" "Filt_ZC_MWave"]);
                     end
                     
                    scale=500;
                    asfMultip=10;
                    offset=10^-5/2;
                    
                    if ExpLabel=="MVCTrials"
                        MeanMAV=mean(MAV_vEMG(end-stim_freq*1:end-stim_freq*0));
                    elseif ExpLabel=="RCCurveTrials"
                        MeanMAV=mean(MAV_vEMG(end-stim_freq*4:end-stim_freq*2));
                    elseif ExpLabel=="RCRampTrials"
                        MeanMAV=mean(MAV_vEMG(end-stim_freq*2:end-stim_freq*0));
                    else
                        MeanMAV=mean(MAV_vEMG(end-stim_freq*7:end-stim_freq*2));
                    end
    
                    asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11]*MeanMAV*asfMultip;
                    [Amp_MAV_vEMG,Clip_MAV_vEMG]=amp_modul(MAV_vEMG,scale,offset,asf); 
                    [Amp_MAV_MWave,Clip_MAV_MWave]=amp_modul(MAV_MWaves,scale,offset,asf); 
                    
                    Label=sprintf("AmpModulFeats%s",FieldSuffix(iField));                    
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                        table(Amp_MAV_vEMG',Clip_MAV_vEMG',Amp_MAV_MWave',Clip_MAV_MWave',...
                        'VariableNames',[ "Amp_MAV_vEMG" "Clip_MAV_vEMG" ...
                        "Amp_MAV_MWave" "Clip_MAV_MWave"]);

                end
            end
        end
    end
    S.(AnaLabel).AnaPar.FieldSuffix=FieldSuffix;
    S.(AnaLabel).AnaPar.FieldNames=FieldNames;
end

% %%Plotting EMG Features 
% clc
% AnaLabel=sprintf('%s_ana',TestFolders(1));
% ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
% TimeRange=[0.1 13];
% TrialNum=[  4];
% cm=lines(3);
% FiltLabel="GS";
% for iTest=1:length(TestFolders)
%     AnaLabel=sprintf('%s_ana',TestFolders(iTest));
%     TestLabel=sprintf('%s_test',TestFolders(iTest));
%     stim_freq=S.(TestLabel).ExpPar.stim_freq;
%     FrameInd=ceil([stim_freq*TimeRange(1): stim_freq*TimeRange(2)]);
% 
%     for iTrial=1:length(TrialNum)
% 
%         TrialLabel=sprintf('Trial_%d',TrialNum(iTrial));
%         AmpMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Amp_MAV_vEMG');
%         AmpClipped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Clip_MAV_vEMG');   
% 
%         figure(1)
%         subplot(length(TestFolders),1,iTest)
%         plot(FrameInd,AmpMAV,'LineWidth',2,'DisplayName',sprintf('AmpMAV, Trial: %d',TrialNum(iTrial)),'Color',cm(1,:))
%         hold on 
%         subplot(length(TestFolders),1,iTest)
%         plot(FrameInd,AmpClipped,'DisplayName',sprintf('AmpClipped, Trial: %d',TrialNum(iTrial)),'Color',cm(2,:))
%         title(TestFolders(iTest))
%     end
% 
% %     title(ttl)
%     legend('Location','NorthWest')
%     grid on
% end


%%DroppedFrames Features

ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
ExpstoAna=ExpLabels([4]);
DroppedFeat=[];
DroppedFeatTrial=[];
for iExp =1:length(ExpstoAna)
    ExpLabel=ExpstoAna(iExp);

    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders{iTest});
        AnaLabel=sprintf("%s_ana",TestFolders{iTest});

        TrialNum=S.(TestLabel).(ExpLabel).NumofTrials;
        BlankLength=S.(AnaLabel).AnaPar.BlankLength;
        FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
        for iTrial=1:TrialNum
            TrialLabel=sprintf('Trial_%d',iTrial);
            DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
            if isempty(DroppedFrames)
                continue;
            end

            if isfield(S.(TestLabel).(ExpLabel), 'TrialsPW')

                TrialsPW=S.(TestLabel).(ExpLabel).TrialsPW;
                dropped=S.(TestLabel).(ExpLabel).dropped;
                if dropped==1 %&& TrialsPW(iTrial) ~= 0
                    for iFilt=1:length(FiltLabels)
                    FiltLabel=FiltLabels(iFilt);
                    
                    % DroppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG;
                    DroppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(:,DroppedFrames);

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

                    Filts=strings(length(Ssc),1)+FiltLabel;
                    Trials=strings(length(Ssc),1)+TrialLabel;
                    Exps=strings(length(Ssc),1)+ExpLabel;
                    Tests=strings(length(Ssc),1)+TestLabel;

                    DroppedFeat=[DroppedFeat;...
                        MAV_vEMG' MedFreq_vEMG' MeanFreq_vEMG' Ssc' Zc' Filts Trials Exps Tests];
                    DroppedFeatTrial=[DroppedFeatTrial; ...
                        MAV_vEMG' MedFreq_vEMG' MeanFreq_vEMG' Ssc' Zc' Filts];

                    S.(AnaLabel).(ExpLabel).DroppedFeat=array2table(DroppedFeat,'VariableNames',...
                        ["MAV" "MedFreq" "MeanFreq" "SSC" "ZC" "Filt" "Trial" "Exp" "Test"]);
                    end
                    S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat=array2table(DroppedFeatTrial,...
                        'VariableNames',["MAV" "MedFreq" "MeanFreq" "SSC" "ZC" "Filt"]);
                    DroppedFeatTrial=[];
                end
            end
        end
    end
end


%%Mean MAV of Stim Only Trials 
clear MAV_Mean

RCMeanTime=[8 10]; % Calculating the means at time [8 10]
FiltLabel="Unfilt";
RampTimeRange=[2 20];
MAV_Mean_tests=[];
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});

    ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    ExpRuns=S.(TestLabel).ExpRuns;
    RampMAV=[];
    MAVVal=[];
    MAVVal_reps=[];

    if ExpRuns(str2double(S.(TestLabel).ExpPar.ExpTable(2,:).('RC')))
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('RC');

        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        MeanFrame=[RCMeanTime(1)*stim_freq RCMeanTime(2)*stim_freq];
        MeanRangeInd=[MeanFrame(1): MeanFrame(2)];

        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials';
        PercentMVCVals=S.(TestLabel).(ExpLabel).PercentMVC';
        PWPoints=S.(TestLabel).(ExpLabel).PWPoints';
        PWofTrials=S.(TestLabel).(ExpLabel).PWTrials';
        IndTrials=[];

        for iPW=1:length(PWPoints)

            Ind_PW(:,iPW)=PWofTrials==PWPoints(iPW);
            IndTrials(:,iPW)=find(Ind_PW(:,iPW)==1);

            for iTrial=1:length(IndTrials(:,iPW))
                TrialLabel=sprintf('Trial_%d',IndTrials(iTrial,iPW));
                MAV_mean_reps(iTrial)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel)...
                    .(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
                MAV_std_reps(iTrial)=std(S.(AnaLabel).(ExpLabel).(TrialLabel)...
                    .(FiltLabel).Feats.('MAV_vEMG')(MeanRangeInd));
                
               MAVVal_reps=[MAVVal_reps; MAV_mean_reps(iTrial) MAV_std_reps(iTrial) PWPoints(iPW)  "RC" TestFolders(iTest) ];

            end
            
           MAVVal=[MAVVal; mean(MAV_mean_reps) std(MAV_mean_reps) PWPoints(iPW)  "RC" TestFolders(iTest) ];

        end
        S.(AnaLabel).(ExpLabel).MAV_Mean_Reps=array2table(MAVVal_reps,...
            'VariableNames',["Mean" "Std" "PW" "Exp" "Test"] );
        MAV_Mean_tests=[MAV_Mean_tests; MAVVal];

    end
    
    if  ExpRuns(double(S.(TestLabel).ExpPar.ExpTable(2,:).('Ramp')))
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Ramp');   
        
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        MVCLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
        OccLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');

        MVC=S.(TestLabel).(MVCLabel).MVC;
%         PercentMVCVals=MVC/100*linspace(min(S.(TestLabel).(OccLabel).VoliMVCVec),max(S.(TestLabel).(OccLabel).VoliMVCVec),length(PWPoints));
        PercentMVCVals=MVC*S.(TestLabel).(OccLabel).VoliMVCVec/100;

        for iTrial=1:NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial);
            TrialFrameLength=length(S.(AnaLabel).(ExpLabel).(TrialLabel).BegofFrames);
            stim_freq=S.(TestLabel).ExpPar.stim_freq;

            FrameRange=stim_freq*RampTimeRange;
            FrameRange(FrameRange>TrialFrameLength)=TrialFrameLength-1;
            FrameRangeInd=FrameRange(1):FrameRange(2);
            
            MAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Filt_MAV_vEMG')(FrameRangeInd);
            PWLevels=S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames(FrameRangeInd);
            Force=S.(AnaLabel).(ExpLabel).(TrialLabel).FiltForceFrames(FrameRangeInd);
                             
            RampMAVTrial=[[],[],[],[],[],[]];
            IndMin=[];
            MeanRampPWVals=[];
            MeanRampMAVVals=[];
            StdRampMAVVals=[];
            MeanRampForceVals=[];
            StdRampForceVals=[];
            for iVal=1:length(PercentMVCVals)
                
                [ValMin(iVal), IndMin(iVal)]=min(abs(PercentMVCVals(iVal)-Force));
                RangeLow=IndMin(iVal)-100;
                RangeHigh=IndMin(iVal)+100;
                RangeInd=RangeLow:RangeHigh;
                RangeInd(RangeInd<1)=[];
                RangeInd(RangeInd>length(MAV))=[];
                
                MeanRampMAVVals(iVal)=mean(MAV(RangeInd));
                StdRampMAVVals(iVal)=std(MAV(RangeInd));
                MeanRampForceVals(iVal)=mean(Force(RangeInd));
                StdRampForceVals(iVal)=std(Force(RangeInd));
                MeanRampPWVals(iVal)=mean(PWLevels(RangeInd));
                RampMAV=[RampMAV; MeanRampMAVVals(iVal) StdRampMAVVals(iVal) MeanRampForceVals(iVal)...
                    StdRampForceVals(iVal) MeanRampPWVals(iVal) IndMin(iVal)+RampTimeRange(1)*stim_freq,...
                    TrialLabel, TestFolders(iTest)];

            end
                        
            S.(AnaLabel).(ExpLabel).(TrialLabel).PercentMVCVals=PercentMVCVals;
            S.(AnaLabel).(ExpLabel).(TrialLabel).ValMin=ValMin;
            S.(AnaLabel).(ExpLabel).(TrialLabel).Frame=IndMin +RampTimeRange(1)*stim_freq;
            S.(AnaLabel).(ExpLabel).(TrialLabel).MeanRampMAV=MeanRampMAVVals;
            S.(AnaLabel).(ExpLabel).(TrialLabel).StdRampMAV=StdRampMAVVals;
            S.(AnaLabel).(ExpLabel).(TrialLabel).MeanRampForceVals=MeanRampForceVals;
            S.(AnaLabel).(ExpLabel).(TrialLabel).StdRampForceVals=StdRampForceVals;     
            S.(AnaLabel).(ExpLabel).(TrialLabel).MeanRampPWVals=MeanRampPWVals;
            
            RampMAVTrial=[RampMAVTrial; MeanRampMAVVals' StdRampMAVVals' MeanRampForceVals'...
                StdRampForceVals' MeanRampPWVals' IndMin'+RampTimeRange(1)*stim_freq];
                
            RampMAVTrialTable=array2table(RampMAVTrial,'VariableNames',["MeanRampMAV","SDRampMAV",...
                "MeanRampForceVals","StdRampForceVals","MeanRampPWVals","Frame"]);
            
            S.(AnaLabel).(ExpLabel).(TrialLabel).RampMAVTable=RampMAVTrialTable;
        end
    
    RampMAVTable=table(RampMAV(:,1),RampMAV(:,2),RampMAV(:,3),RampMAV(:,4),RampMAV(:,5),RampMAV(:,6),RampMAV(:,7),RampMAV(:,8),...
    'VariableNames',["MeanRampMAV","SDRampMAV","MeanRampForceVals","StdRampForceVals","MeanRampPWVals","Frame", "Trial", "Test"]);
    S.(AnaLabel).(ExpLabel).RampMAVTable=RampMAVTable;
    
    MAVVal=[MAVVal; MeanRampMAVVals' StdRampMAVVals' (MeanRampPWVals)' strings(length(MeanRampMAVVals),1)+"Ramp"...
        strings(length(MeanRampMAVVals),1)+TestFolders(iTest) ]; % Use the last trial
    
    MAV_Mean_tests=[MAV_Mean_tests; MAVVal];
    
    clear MeanRampMAVVals StdRampMAVVals MeanRampPWVals
    
    MAV_Mean=array2table(MAVVal,'VariableNames',["Mean" "Std" "PW" "Exp" "Test"] );
    S.(AnaLabel).(ExpLabel).MAV_Mean=MAV_Mean;
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('RC');
    S.(AnaLabel).(ExpLabel).MAV_Mean=MAV_Mean;
    end 
end

% Update this to more general
MAV_Mean_tests=array2table(MAV_Mean_tests,'VariableNames',["Mean" "Std" "PW" "Exp" "Test"] );

%%Mean MAV of Occ Trials 
clc
clear MAVDropped
exp_lbl='Occ';
% VoliLevels=[10 20 30 40];
% StimLevels=[ 10 12 15 0];
TrialStats=[];

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(exp_lbl);
    VoliLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
    StimLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    
    MeanTime=S.(TestLabel).(ExpLabel).StimConstantRange;
    MeanFrame=[MeanTime(1)*stim_freq MeanTime(2)*stim_freq];
    MeanRangeInd=[MeanFrame(1):MeanFrame(2)];
    c=1;
    clear MAV_Mean MAV_Mean_Reps TrialsInd Amp_Modul_Mean  Amp_Modul_Mean_Reps Mean_MAVDropped Std_MAVDropped RepTrialsInd 
    TrialStats=[];
    Mean_TrialStats=[];
    for iVoli=1:length(VoliLevels)
        VoliLevel=VoliLevels(iVoli);
        FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;

        for iStim=1:length(StimLevels)
            StimLevel=StimLevels(iStim);

            TrialsInd(:,c)= find_trialnum(VoliLevel, StimLevel, RepTableMat);
            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels(iFilt);

                for iTrial=1:length(TrialsInd(:,c))
                    TrialLabel=sprintf('Trial_%d',TrialsInd(iTrial,c));

                    if StimLevel==0
                        Mean_MAVDropped(TrialsInd(iTrial,c),1)=NaN;
                        Std_MAVDropped(TrialsInd(iTrial,c),1)=NaN;
                    
                    else 
                        DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
                        FiltInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat.('Filt')==FiltLabel;
                        DropFrames=double(S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(FiltInd,:).('MAV'));
    
                        Mean_MAVDropped(TrialsInd(iTrial,c),1)=mean(DropFrames(MeanFrame(1)<=DroppedFrames & DroppedFrames<=MeanFrame(2))); 
                        Std_MAVDropped(TrialsInd(iTrial,c),1)=std(DropFrames(MeanFrame(1)<=DroppedFrames & DroppedFrames<=MeanFrame(2)));    
                    end
                
                    MAV_Mean_reps(TrialsInd(iTrial,c),1)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(MeanRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),2)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(MeanRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),3)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(MeanRangeInd,:).('Amp_MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),4)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(MeanRangeInd,:).('Amp_MAV_vEMG'));                    
                    MAV_Mean_reps(TrialsInd(iTrial,c),5)=Mean_MAVDropped(TrialsInd(iTrial,c),1);
                    MAV_Mean_reps(TrialsInd(iTrial,c),6)=Std_MAVDropped(TrialsInd(iTrial,c),1);
                    MAV_Mean_reps(TrialsInd(iTrial,c),7:9)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);

                    Amp_Modul_Mean(TrialsInd(iTrial,c),3:5)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);

                    TrialStats=[TrialStats; MAV_Mean_reps(TrialsInd(iTrial,c),:) FiltLabel TrialsInd(iTrial,c) ExpLabel TestLabel  ];                               
                end

                MAV_Mean(c,:) = [ c mean(MAV_Mean_reps(TrialsInd(:,c),1))  std(MAV_Mean_reps(TrialsInd(:,c),1)) mean(MAV_Mean_reps(TrialsInd(:,c),3))...
                    std(MAV_Mean_reps(TrialsInd(:,c),3)) mean(Mean_MAVDropped(TrialsInd(:,c),1)) std(Mean_MAVDropped(TrialsInd(:,c),1))  ...
                    length(TrialsInd(:,c)) RepTableMat(TrialsInd(iTrial,c),[1 4 5 ])];

                Mean_TrialStats=[Mean_TrialStats; MAV_Mean(c,:) FiltLabel TrialsInd(iTrial,c) ExpLabel TestLabel];
            end
    
            c=c+1;
        end
    end
    
    MAV_Mean_reps_table=table(double(TrialStats(:,1)),double(TrialStats(:,2)),double(TrialStats(:,3)),double(TrialStats(:,4)),double(TrialStats(:,5)),...
        double(TrialStats(:,6)),double(TrialStats(:,7)),double(TrialStats(:,8)),double(TrialStats(:,9)),TrialStats(:,10),TrialStats(:,11),TrialStats(:,12),...
        TrialStats(:,13),'VariableNames',["MAV_Mean" "MAV_Std" "Amp_Mean" "Amp_Std" "Dropped_Mean" "Dropped_Std" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "Exp" "Test"]);
    
    MAV_Mean_table=table(double(Mean_TrialStats(:,1)),double(Mean_TrialStats(:,2)),double(Mean_TrialStats(:,3)),double(Mean_TrialStats(:,4)),double(Mean_TrialStats(:,5)),...
        double(Mean_TrialStats(:,6)),double(Mean_TrialStats(:,7)),double(Mean_TrialStats(:,8)),double(Mean_TrialStats(:,9)),double(Mean_TrialStats(:,10)),double(Mean_TrialStats(:,11)),...
        Mean_TrialStats(:,12),Mean_TrialStats(:,13),Mean_TrialStats(:,14),Mean_TrialStats(:,15),'VariableNames',["Unique_Trial" "MAV_Mean_Reps" "MAV_Std_Reps"...
        "Amp_Mean_Reps" "Amp_Std_Reps" "Dropped_Mean_Reps" "Dropped_Std_Reps" "NumofReps" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "Exp" "Test"]);
    
    S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table=MAV_Mean_reps_table;
    S.(AnaLabel).(ExpLabel).MAV_Mean_table=MAV_Mean_table;
    
end

%% Plotting
% 
% cm = lines(length(TestFolders));
% sMVC=0;
% FiltLabel="Unfilt";
% Extrapolate=100;
% for iTest=1:length(TestFolders)
%     AnaLabel=sprintf("%s_ana",TestFolders{iTest});
%     ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
%     FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');
% 
%     sMVCzero_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC');
%     MAVMean_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean');
%     vMVCVal_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC');
% 
%     % sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
%     % MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps');
%     % vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
%     % Ind=(sMVCzero==sMVC);
% 
%     Amp_Modul_Mean=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean');
% 
%     Ind_Reps=(sMVCzero_Reps==sMVC);
%     p1=polyfit(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),1);
% 
%     p2=polyfit(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),1);
% 
%     figure(1)
%     subplot(2,1,1)
%     plot(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
%     legend
%     hold on
%     plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p1,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
%     xlabel('% Voli. MVC')
%     ylabel('Mean MAV')
%     % xlim([0 50])
% 
%     subplot(2,1,2)
%     plot(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
%     legend
%     hold on
%     plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p2,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
%     ylabel('% Amp modul')
%     xlabel('% Voli. MVC')
%     % xlim([0 50])
% end

%% Normalize the Features over zero stim trials
FeatLabels=string(S.(AnaLabel).AnaPar.FeatLabels);
% VoliMVCLevels=[10 20 30 40 ];
% StimMVCLevels=[0 10 20 30 ];
FiltLabel="Unfilt";

for iFeat=1:1
    FeatLabel=FeatLabels(iFeat);
    
    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders{iTest});
        AnaLabel=sprintf("%s_ana",TestFolders{iTest});

        FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_table.('Filt');

        RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
        VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
        StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;

        for iVoli=1:length(VoliMVCLevels)
            sMVC_Vec=(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('sMVC'));
            vMVC_Vec=(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('vMVC'));

            Voli_NormCoeff=(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('MAV_Mean_Reps')...
                (sMVC_Vec==0 & vMVC_Vec==VoliMVCLevels(iVoli)));

            for iStim=1:length(StimMVCLevels)
                sMVC=StimMVCLevels(iStim);
                vMVC=VoliMVCLevels(iVoli);
                IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

                for iRep=1:length(IndTrials)
                    TrialLabel=sprintf('Trial_%d',IndTrials(iRep));
                    FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;

                    for iFilt=1:length(S.(AnaLabel).AnaPar.FiltLabels)
                        FiltLabel=S.(AnaLabel).AnaPar.FiltLabels{iFilt};
                        vEMGLabel=sprintf('%s_vEMG',FeatLabel);
                        MWaveLabel=sprintf('%s_MWave',FeatLabel);
                        
                        FeatvEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel);
                        FeatMWave=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(MWaveLabel);

                        NormvEMGFeat=FeatvEMG/Voli_NormCoeff;
                        NormMWaveFeat=FeatMWave/Voli_NormCoeff;

                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats...
                            .(sprintf('Norm_%s',vEMGLabel))=NormvEMGFeat;
                        
                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats...
                            .(sprintf('Norm_%s',MWaveLabel))=NormMWaveFeat;

                    end
                end
            end
        end
    end
end

%% Theoretical and Actual MVC MAV 

NumofUpdate=4; % table mvc_table is updated 4 times inside loop below
NumofVariables=4;
mvc_table=strings(NumofUpdate*length(TestFolders),NumofVariables);
MVC_Percent=100;
sMVC_Ref=0;
AvgTimeMVC=2;
FiltLabel="Unfilt";
cm = lines(length(TestFolders));

for iTest=1:length(TestFolders)
    
    % MAV_MAX = avg MAV at the last 2 seconds
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    stim_freq=S.(TestLabel).ExpPar.stim_freq;

    AvgInd=stim_freq*AvgTimeMVC-10;

    TrialsMAV=zeros(1,S.(TestLabel).(ExpLabel).NumofTrials);
    TrialsAmp_MAV=zeros(1,S.(TestLabel).(ExpLabel).NumofTrials);
    for iTrial=1:S.(TestLabel).(ExpLabel).NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial );

        TrialsMAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
            .Feats(end-AvgInd:end,:).('MAV_vEMG'));
        
        TrialsAmp_MAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
            .AmpModulFeats(end-AvgInd:end,:).('Amp_MAV_vEMG'));
    end
    
    MAV_max=max(TrialsMAV);
    S.(AnaLabel).(ExpLabel).MAV_MAX=MAV_max;
    
    AmpModul_MAV_max=max(TrialsAmp_MAV);
    S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX=AmpModul_MAV_max;
    
    mvc_table(iTest*NumofUpdate-3,:)=[MAV_max "MAV" false TestFolders(iTest) ];        
    mvc_table(iTest*NumofUpdate-2,:)=[AmpModul_MAV_max "Amp" false TestFolders(iTest) ];
    
    % Theoretical MAV_MAX
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    sMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC'));
    vMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC'));
    
    MAVMean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean'));
    Amp_Modul_Mean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean'));
    
    Ind=(sMVCLevs==sMVC_Ref);
    
    poly1=polyfit(vMVCLevs(Ind),MAVMean(Ind),1);
    poly2=polyfit(vMVCLevs(Ind),Amp_Modul_Mean(Ind),1);
    
    MAV_max_theo=polyval(poly1,MVC_Percent);
    AmpModul_max_theo=polyval(poly2,MVC_Percent);
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    S.(AnaLabel).(ExpLabel).MAV_MAX_theo=MAV_max_theo;
    S.(AnaLabel).(ExpLabel).AmpModul_max_theo=AmpModul_max_theo;
    
    mvc_table(iTest*NumofUpdate-1,:)=[ MAV_max_theo "MAV" true TestFolders(iTest) ];
    mvc_table(iTest*NumofUpdate,:)=[ AmpModul_max_theo "Amp" true TestFolders(iTest) ];
    
    S.(AnaLabel).(ExpLabel).MVCTable= array2table(mvc_table,'VariableNames',["MVC_MAV" "Type" "Theo" "Test"]);


%%-----> Plotting the MVC_ MAV levels

    % figure(1) 
    % subplot(2,1,1)
    % semilogy(vMVCLevs(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    % hold on 
    % semilogy([unique(vMVCLevs(Ind)); 100],polyval(poly1,[unique(vMVCLevs(Ind)); 100]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest));
    % semilogy(100, MAV_max, '*','Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    % legend
    % xlabel('% Voli. MVC')
    % ylabel('Mean MAV')
    % % xlim([0 50])
    
end

writetable(S.(AnaLabel).(ExpLabel).MVCTable,'mvc_table.csv')




%% Occlusion Analysis

% AvgOcclusionTests=[  "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
close all
% TauTests=["jan7" "jan11" "jan12"]; %% Tests for time constant calculation
% TauTests=["feb29_24" "mar18_24"  "mar20_24"]
TestIndForce=Tau_stats_table.('Test')==TauTestsForce;
TestIndHand=Tau_stats_table.('Test')==TauTestsHand;

taus=Tau_stats_table(any(TestIndForce'),:).('Mean');
MeanTaus=mean(taus); % average time constant to be used 

OccRefs = [3*MeanTaus 4*MeanTaus 5*MeanTaus]; % referance time for occlusion to be calculated
TauLabels=["3*tau" "4*tau" "5*tau"];
OccType=["Occ_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc"];
NoiseAvgTime=[2 3];

Occ={};
for iTau=1:length(OccRefs)
    OccRef=OccRefs(iTau);

    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders(iTest));
        AnaLabel=sprintf("%s_ana",TestFolders(iTest));
        FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');

        EffortType=S.(TestLabel).ExpPar.EffortType;
        DroppedFiltLabel=DroppedFiltLabels(iTest);
        NoStimFiltLabel=NoStimFiltLabels(iTest);
        stim_freq=S.(TestLabel).ExpPar.stim_freq;

        if MAV_MAXMethods(iTest)=="Fitted"
            MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;           %% MAV_MAX replacement 
            AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;
        else
            MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX;
            AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;
        end

        MVC=S.(TestLabel).(ExpLabel).MVC;

        AvgInd=round(stim_freq*NoiseAvgTime(1):stim_freq*(NoiseAvgTime(2)));
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        for iTrial=1:S.(TestLabel).(ExpLabel).NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial );

            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels(iFilt);
                MAV_Noise= mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(AvgInd,:).('Filt_MAV_vEMG'));
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MAV_Noise=MAV_Noise;

                AmpModul_MAV_Noise= mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(AvgInd,:).('Amp_MAV_vEMG'));
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModul_MAV_Noise=AmpModul_MAV_Noise;
            end
        end
        
        OccRefMargin= 0.2/stim_freq ; % in secs ------------> check if this is correct
        PreOffMargin=10/stim_freq; % Avg period before stim turning off 

        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat(:,1:7),...  %%-- Update this repmattable issue
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);
%         S.(TestLabel).(ExpLabel).RepTableMat=RepMatTable;       


        for iTrial=1:NumofTrials
            TrialLabel=sprintf("Trial_%d", iTrial);
            
            Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
            PostOffInd= TurnOffTime+OccRef-OccRefMargin<=Time & Time <=TurnOffTime+OccRef+OccRefMargin;
            PreOffInd=TurnOffTime-PreOffMargin<=Time & Time <=TurnOffTime;
            
            EffortMeasure=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortType);

            [~,TurnOffInd]=min(abs(Time-TurnOffTime));
            PostForceLevel= mean(EffortMeasure( PostOffInd));
            PreForceLevel= mean(EffortMeasure( PreOffInd));
            PW_level=RepMatTable(iTrial,:).('PW');
            F_v=RepMatTable(iTrial,:).('Voli_Force');
            F_stim=RepMatTable(iTrial,:).('Stim_Force');
            F_target=RepMatTable(iTrial,:).('Target_Level');
            
            F_vprime=PostForceLevel;
            F_sprime=PreForceLevel-PostForceLevel;
            F_occ=F_vprime-F_v;
            
            F_occ_mvc=F_occ/MVC*100;
            F_vprime_mvc=F_vprime/MVC*100;
            F_sprime_mvc=F_sprime/MVC*100;
            F_v_mvc=F_v/MVC*100;
            F_target_mvc=F_target/MVC*100;
            
            PostOffFrameInd=round((TurnOffTime+OccRef-OccRefMargin)*stim_freq:( TurnOffTime+OccRef+OccRefMargin)*stim_freq);
            PreOffFrameInd=round((TurnOffTime-PreOffMargin)*stim_freq:( TurnOffTime)*stim_freq);
            % FiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('Filt')==FiltLabel;
            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('Filt')==NoStimFiltLabel;
            FiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==DroppedFiltLabel;
            TrialNums_reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Trial');

            MVC_Voli=RepMatTable(iTrial,:).('MVC_Voli');
            MVC_Stim=RepMatTable(iTrial,:).('MVC_Stim');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
            MAV_v=S.(AnaLabel).(ExpLabel).MAV_Mean_table(NoStimFiltInd & MVC_Voli==vMVC & sMVC==0,:).('MAV_Mean_Reps'); % Are those means before or after the turn off ????
            AmpModul_MAV_v=S.(AnaLabel).(ExpLabel).MAV_Mean_table(NoStimFiltInd & MVC_Voli==vMVC & sMVC==0,:).('Amp_Mean_Reps');

            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==NoStimFiltLabel;
            DroppedMean_temp=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Dropped_Mean');

            DroppedMAV=DroppedMean_temp(iTrial==double(TrialNums_reps),:);
            vMVC_temp=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd,:).('vMVC');
            vMVC_iTrial=vMVC_temp(iTrial==double(TrialNums_reps),:);

            vMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('vMVC');
            sMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('sMVC');
            Mean_NoStimMAV=mean(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd & vMVCTrials==vMVC_iTrial & sMVCTrials==0,:).('MAV_Mean'));
            
            % 1- E_o = E_v-E_d       (Dropped) 
            % 2- E_o = (E_v-E_d)/mvc*100       (Dropped Effort) 
            % 3- E_o = (E_vprime-E_d)/mvc*100   (Hybrid Effort)

            Dropped_occ= DroppedMAV-Mean_NoStimMAV;
            Dropped_occ_mvc= Dropped_occ/MAV_max*100;
            Hybrid_occ_mvc= F_vprime_mvc-Dropped_occ_mvc;
            
            Occ=[ Occ;  F_occ  F_vprime F_sprime F_v Dropped_occ F_target F_occ_mvc F_vprime_mvc...
                F_sprime_mvc F_v_mvc Dropped_occ_mvc Hybrid_occ_mvc DroppedMAV Mean_NoStimMAV MAV_max...
                F_target_mvc TestFolders(iTest) EffortType TauLabels(iTau) "EffortMea" "Unfilt" ...
                table2array(RepMatTable(iTrial,:) ) iTrial ];

                %%s MAV based occlusion estimation using stim-turned off region 
            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels{iFilt};

                MAV_target=F_target; %%% ---- Fix This 

                PreOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(PreOffFrameInd,:).('MAV_vEMG'));
                PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(PostOffFrameInd,:).('MAV_vEMG'));
                MAV_vprime=PostOffMAV;
                MAV_sprime=PreOffMAV-PostOffMAV;

                MAV_occ=MAV_vprime-MAV_v;
                MAV_occ_mvc=MAV_occ/MAV_max*100;
                MAV_vprime_mvc=MAV_vprime/MAV_max*100;
                MAV_sprime_mvc=MAV_sprime/MAV_max*100;

                MAV_v_mvc=MAV_v/MAV_max*100;
                MAV_target_mvc=MAV_target/MAV_max*100;

                Occ=[ Occ;  MAV_occ  MAV_vprime MAV_sprime MAV_v Dropped_occ MAV_target MAV_occ_mvc ...
                    MAV_vprime_mvc MAV_sprime_mvc MAV_v_mvc Dropped_occ_mvc Hybrid_occ_mvc  DroppedMAV Mean_NoStimMAV MAV_max MAV_target_mvc...
                    TestFolders(iTest) EffortType TauLabels(iTau) "MAV" FiltLabel table2array(RepMatTable(iTrial,:)) iTrial];

                AmpModul_MAV_target=F_target; %%% ---- Fix This 

                AmpModul_PreOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                    .AmpModulFeats(PreOffFrameInd,:).('Amp_MAV_vEMG'));
                AmpModul_PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                    .AmpModulFeats(PostOffFrameInd,:).('Amp_MAV_vEMG'));

                AmpModul_MAV_vprime=AmpModul_PostOffMAV;
                AmpModul_MAV_sprime=AmpModul_PreOffMAV-AmpModul_PostOffMAV;

                AmpModul_MAV_occ=AmpModul_MAV_vprime-AmpModul_MAV_v;
                AmpModul_MAV_occ_mvc=AmpModul_MAV_occ/AmpModul_MAV_max*100;
                AmpModul_MAV_vprime_mvc=AmpModul_MAV_vprime/AmpModul_MAV_max*100;
                AmpModul_MAV_sprime_mvc=AmpModul_MAV_sprime/AmpModul_MAV_max*100;

                AmpModul_MAV_v_mvc=AmpModul_MAV_v/AmpModul_MAV_max*100;
                AmpModul_MAV_target_mvc=AmpModul_MAV_target/AmpModul_MAV_max*100;

                Occ=[ Occ; AmpModul_MAV_occ AmpModul_MAV_vprime AmpModul_MAV_sprime AmpModul_MAV_v ...
                    Dropped_occ AmpModul_MAV_target AmpModul_MAV_occ_mvc AmpModul_MAV_vprime_mvc...
                    AmpModul_MAV_sprime_mvc AmpModul_MAV_v_mvc Dropped_occ_mvc Hybrid_occ_mvc  DroppedMAV Mean_NoStimMAV MAV_max  ...
                    AmpModul_MAV_target_mvc TestFolders(iTest) EffortType TauLabels(iTau) "Amp_Modul"...
                    FiltLabel table2array(RepMatTable(iTrial,:)) iTrial];
            end
        end

        OccTest=Occ(Occ(:,14)==TestFolders(iTest),:); % 14 is where Test variable is located
        OccTest=array2table(OccTest,'VariableNames',[ "Occ" "vprime" "sprime" "v" "Occ_Dropped" "Target" "Occ_mvc" "vprime_mvc"...
            "sprime_mvc" "v_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc" "DroppedMAV" "NoStimMAV" "MAV_max" "Target_mvc" "Test" ...
            "EffortType" "Tau" "Feat" "Filt" "Target_Level" "Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW" "Done" "Trial"]);
        
        S.(AnaLabel).(ExpLabel).OccTest=OccTest;
        clear OccTest
    end
end

OccTable=array2table(Occ,'VariableNames', [ "Occ" "vprime" "sprime" "v" "Occ_Dropped" "Target" "Occ_mvc" "vprime_mvc"...
            "sprime_mvc" "v_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc" "DroppedMAV" "NoStimMAV" "MAV_max" "Target_mvc" "Test" ...
            "EffortType" "Tau" "Feat" "Filt" "Target_Level" "Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW" "Done" "Trial"]);

writetable( OccTable,'occlusion_v7.csv')


%% Occlusion Fitting  
%linear modeling for individual occ predictions 
    % Individual occlusion estimation
clear lm_table
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20"];
lm_table=table();
CoefMat=[];
CoefMat2=[];
CoefMat3=[];

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    EffortTypeInd(iTest)=S.(TestLabel).ExpPar.EffortType;
end

EffortTypes=unique(EffortTypeInd);
for iEffortType=1:length(EffortTypes)
    EffortLabel=EffortTypes(iEffortType);
    EffortTestFolders=TestFolders(EffortTypeInd==EffortLabel);
    Coefs=[];
    for iType=1:length(OccType)
       lm_table=table();
    
        for iTest=1:length(EffortTestFolders)
            TestLabel=sprintf("%s_test",EffortTestFolders(iTest));
    
            EffortType=S.(TestLabel).ExpPar.EffortType;
            RowInd=OccTable.('EffortType')==EffortType & OccTable.('Feat')=="EffortMea" & OccTable.('Filt')=="Unfilt"...
                & OccTable.('Tau')=="3*tau" & OccTable.('MVC_Stim')~="0";
            
    % 1- Linear fitting for Individualized results
            lm_table.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
            lm_table.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
            lm_table.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
            lm_table.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
            lm_table.Test=categorical((OccTable(RowInd,:).('Test')));
            
            lm_table.Test=reordercats(lm_table.Test,EffortTestFolders);
            mdl = fitlm(lm_table,'Occ~VoliMVC+StimMVC+Test');
            CoefNum=table2array(mdl.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'))';
            CoefCat=table2array(mdl.Coefficients([1,4:end],'Estimate'));
            pVals=mdl.Coefficients.('pValue');
            pValInd=[1 4:4+length(EffortTestFolders)-2]; % 1: p val for first test, 2,3: numerical variables, 4-:rest of the tests  
            TestInd=lm_table.Test==EffortTestFolders(iTest);
    
            if iTest==1
                Coefs(iTest,:)=[CoefCat(1) CoefNum];
            else
                Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum];
            end
            
            Effort_o=[];
            Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC(TestInd)+Coefs(iTest,3)*lm_table.StimMVC(TestInd))];
            
            CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) pVals(2) pVals(3) pVals(pValInd(iTest)) EffortTestFolders(iTest) EffortType...
                OccType(iType) boolean(1) boolean(0)], length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd)...
                lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]'];
            
            CoefMat3=[ CoefMat3; Coefs(iTest,:) pVals(pValInd(iTest)) pVals(2) pVals(3) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(0) ];
            
        end
        
        CoefMat=[ CoefMat; Coefs EffortTestFolders' strings(length(EffortTestFolders),1)+EffortType strings(length(EffortTestFolders),1)+OccType(iType) boolean(ones(length(EffortTestFolders),1))  boolean(zeros(length(EffortTestFolders),1))];
        mdlLog = fitlm(lm_table,'LogOcc~VoliMVC+StimMVC+Test');
        CoefNum=table2array(mdlLog.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
        CoefCat=table2array(mdlLog.Coefficients([1,4:end],'Estimate'));
        pVals=mdlLog.Coefficients.('pValue');
        pValInd=[1 4:4+length(EffortTestFolders)-2]; % 1: p val for first test, 2,3: numerical variables, 4-:rest of the tests  
    %     Effort_o=[];
        for iTest=1:length(EffortTestFolders)
            TestInd=lm_table.Test==EffortTestFolders(iTest);
    
            if iTest==1
                Coefs(iTest,:)=[CoefCat(1) CoefNum'];
            else
                Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum'];
            end
            
            Effort_o=[];
            Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC(TestInd)+Coefs(iTest,3)*lm_table.StimMVC(TestInd))];
            
            CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) pVals(2) pVals(2) pVals(pValInd(iTest)) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(1)],...
                length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd)...
                lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
            
            CoefMat3=[ CoefMat3; Coefs(iTest,:) pVals(pValInd(iTest))  pVals(2) pVals(3) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(1) ];
        end
        
        CoefMat=[ CoefMat; Coefs  EffortTestFolders' strings(length(EffortTestFolders),1)+EffortType strings(length(EffortTestFolders),1)+OccType(iType) boolean(ones(length(EffortTestFolders),1)) boolean(ones(length(EffortTestFolders),1))];
        
    %----- 2- Linear fitting for Averaged results
        lm_table.Test=string(lm_table.Test);
        TestInd=any(lm_table.('Test') == AvgOcclusionTests(iEffortType,:),2);
        lm_table_gen=lm_table(TestInd,:);
    
        mdl_gen = fitlm(lm_table_gen,'Occ~VoliMVC+StimMVC');
        CoefNum_gen=table2array(mdl_gen.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
        CoefCat_gen=table2array(mdl_gen.Coefficients([1,4:end],'Estimate'));
        pVals=mdl_gen.Coefficients.('pValue');
    
    %     Effort_o=[];
    %     Effort_o=[ Effort_o; CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC+CoefNum_gen(2)*lm_table.StimMVC];
    
        CoefMat=[ CoefMat;  CoefCat_gen CoefNum_gen' "Avg" EffortType OccType(iType) boolean(0) boolean(0)];
    
        for iTest=1:length(EffortTestFolders)  
            TestInd=lm_table.Test==EffortTestFolders(iTest);
            Effort_o=[];
            Effort_o=CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC(TestInd)+CoefNum_gen(2)*lm_table.StimMVC(TestInd);
    
            CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(0)],length(lm_table.Occ(TestInd)),1)...
                Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
            
            CoefMat3=[ CoefMat3; CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(0) ];
        end
        
        mdl_genLog = fitlm(lm_table_gen,'LogOcc~VoliMVC+StimMVC');
        CoefNum_gen=mdl_genLog.Coefficients({'VoliMVC', 'StimMVC'},:).('Estimate');
        CoefCat_gen=mdl_genLog.Coefficients([1,4:end],:).('Estimate');
        pVals=mdl_genLog.Coefficients.('pValue');
    
    %     Effort_o=[];
    %     Effort_o=[ Effort_o; CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC+CoefNum_gen(2)*lm_table.StimMVC];
        
        CoefMat=[ CoefMat;  CoefCat_gen CoefNum_gen' "Avg" EffortType OccType(iType) boolean(0) boolean(1)];
        
        for iTest=1:length(EffortTestFolders)  
            TestInd=lm_table.Test==EffortTestFolders(iTest);
            Effort_o=[];
            Effort_o=CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC(TestInd)+CoefNum_gen(2)*lm_table.StimMVC(TestInd);
    
            CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(1)],length(lm_table.Occ(TestInd)),1)...
                Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
            
            CoefMat3=[ CoefMat3; CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(1)];
        end
    end
end


OccCoef=array2table(CoefMat,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "Test" "EffortType" "Type" "Indiv" "Log"]);

OccCoef2=array2table(CoefMat2,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "pVal1" "pVal2" "pVal3" "Test" "EffortType" "Type" "Indiv" "Log"...
    "Effort_o" "Occ" "LogOcc" "MVC_Voli" "MVC_Stim" "Trial"]);

OccCoef3=array2table(CoefMat3,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "pVal1" "pVal2" "pVal3" "Test" "EffortType" "Type" "Indiv" "Log"]);
writetable( OccCoef3,'occ_coef3.csv')


%% Effort Estimation
%Average Occlusion estimation
% OccCoef=[-1.7 -0.109 0.919];  % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fs_primeCoef=[2.117 0.061 0.118 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fv_primeCoef=[-1.361 0.891 0.848 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]
%
% Implementing the correction for occlusion

OccKickOffLevel=0.3;  % Percent effort
ErrMat=[];
LogModel= "false" ;
FeattoAna="Filt_MAV_vEMG";
        % RampTime=[5 10];
        % ConstTime=[10 15];
        % TurnOffTime=[11 12];
clear Effort
for iModel=1:length(LogModel)
    for iTest=1:length(TestFolders)
        AnaLabel=sprintf('%s_ana',TestFolders(iTest));
        TestLabel=sprintf("%s_test",TestFolders(iTest));
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        EffortType=S.(TestLabel).ExpPar.EffortType;
    
        StimRange=S.(TestLabel).(ExpLabel).StimRange;
        StimConstantRange=S.(TestLabel).(ExpLabel).StimConstantRange;
        StimProfile=S.(TestLabel).(ExpLabel).StimProfile;
        RampRange=[StimRange(1) StimConstantRange(1)];

        RampInd=stim_freq*RampRange(1):stim_freq*RampRange(2);
        ConstInd=stim_freq*StimConstantRange(1):stim_freq*StimConstantRange(2);
        FrameInd=[ RampInd ConstInd];

        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
        if MAV_MAXMethods(iTest)=="Fitted"
            MAV_MAX=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;
            AmpModul_MAV_MAX=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;
        else
            MAV_MAX=S.(AnaLabel).(ExpLabel).MAV_MAX;
            AmpModul_MAV_MAX=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;        
        end
        F_MAX=S.(TestLabel).(ExpLabel).MVC;

        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
        OccTest=S.(AnaLabel).(ExpLabel).OccTest; 
        ExpRuns=S.(TestLabel).ExpRuns;
        if ExpRuns(double(S.(TestLabel).ExpPar.ExpTable(2,:).('RC')))
            RCLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('RC');
            RCVar=S.(TestLabel).(RCLabel).RCVar;   
        elseif ExpRuns(double(S.(TestLabel).ExpPar.ExpTable(2,:).('Ramp')))
            Label=S.(TestLabel).ExpPar.ExpTable(1,:).('Ramp');
            RCVar=S.(TestLabel).(Label).RCVar;   
        end
        % What if there is no RCVAR at all ?
        MVCLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
        MVC=S.(TestLabel).(MVCLabel).MVC; 

        RepTableMat=array2table(S.(TestLabel).(ExpLabel).RepTableMat(:,1:7),...  %- Fix this repmattable 
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);

        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        for iTrial=1:NumofTrials 
            TrialLabel=sprintf("Trial_%d",iTrial);

            DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
            PWofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames;
            StimMVC=S.(TestLabel).(ExpLabel).RepTableMat(iTrial,5);
            StimMVCofFrames=gompertz(RCVar,PWofFrames)/MVC*100+0.00001;

            KeepInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd(1:end-1);
            EffortTypeLabel=sprintf("%sFrames",EffortType);

            Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(EffortTypeLabel)(:,KeepInd))';
            Effort_f=Force/F_MAX*100;
            TargetFrames=mean(S.(AnaLabel).(ExpLabel).TargetFrames(:,KeepInd))';
            Target_mvc=RepTableMat(iTrial,:).('Target_Level')/F_MAX*TargetFrames;
            VoliEffort=RepTableMat(iTrial,:).('MVC_Voli')*TargetFrames;
            StimEffort=RepTableMat(iTrial,:).('MVC_Stim')*TargetFrames;

            IndTrials=find_trialnum(RepTableMat(iTrial,:).('MVC_Voli'), ...
                RepTableMat(iTrial,:).('MVC_Stim'),S.(TestLabel).(ExpLabel).RepTableMat);

            RowIndTemp=OccTable.('EffortType')==EffortType & OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')== "EffortMea" & ...
            OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau" & OccTable.('Trial')==string(IndTrials');

            [~,c]=size(RowIndTemp);
            RowInd1=RowIndTemp(:,1);
            for i=1:c-1
                RowInd1=RowInd1 | RowIndTemp(:,i+1);
            end
            %sum(RowInd)

            Fv_prime=str2double(OccTable(RowInd1,:).('vprime_mvc'));

            Effort_Fv_prime=TargetFrames*mean(Fv_prime);

            for iType=1:length(OccType)
                TypeLabel=OccType(iType);

                RowInd=OccCoef3.('Type') == OccType(iType) & OccCoef3.('Test') == TestFolders(iTest) & OccCoef3.('Log') == LogModel(iModel)...
                    & OccCoef3.('Indiv') == "false";
                                
                CoefAvg=str2double([ OccCoef3(RowInd,:).('Coeff1') OccCoef3(RowInd,:).('Coeff2') OccCoef3(RowInd,:).('Coeff3') ]);
                Effort_o=(CoefAvg(1)+CoefAvg(2).*VoliEffort+CoefAvg(3).*StimEffort);

                RowInd=OccCoef3.('Type') == OccType(iType) & OccCoef3.('Test') == TestFolders(iTest) & OccCoef3.('Log') == LogModel(iModel)...
                    & OccCoef3.('Indiv') == "true";

                CoefIndv=str2double([ OccCoef3(RowInd,:).('Coeff1') OccCoef3(RowInd,:).('Coeff2') OccCoef3(RowInd,:).('Coeff3') ]);
                Effort_o_Indv=(CoefIndv(1)+CoefIndv(2).*VoliEffort+CoefIndv(3).*StimEffort);
                
                if LogModel(iModel)=="true"
                    Effort_o=exp(Effort_o);
                    Effort_o_Indv=exp(Effort_o_Indv);
                end
                
                for iFilt=1:length(FiltLabels) 
                    FiltLabel=FiltLabels(iFilt);

                    Filt_Feats=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Filt_MAV_vEMG');
                    MAV_Noise=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MAV_Noise;
                    Filt_Feats=Filt_Feats(1:end-1);                                              %% ---_Fix This 
                        % Need error against dropped frames 
                    % 2-MAV
                        Effort_e_MAV=(Filt_Feats-MAV_Noise)/(MAV_MAX-MAV_Noise)*100;
%                     Effort_e_MAV=Filt_Feats/MAV_MAX;
                    Effort_MAV=Effort_e_MAV+Effort_o;
                    Effort_MAV_Indv=Effort_e_MAV+Effort_o_Indv;

                    Effort_e_MAV_err=Effort_e_MAV-Effort_Fv_prime;
                    Effort_MAV_err=Effort_MAV-Effort_Fv_prime;
                    Effort_MAV_Indv_err=Effort_MAV_Indv-Effort_Fv_prime;

    %                     stdError_MAV=Effort_err_MAV./Fv_prime_frames;
    %                     MeanStdError_MAV=mean(stdError_MAV(2+FrameInd));

    %                     % 3-Amp Modul
                    AmpModul_Feats=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Amp_MAV_vEMG');
                    AmpModul_MAV_Noise=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModul_MAV_Noise;
                    AmpModul_Feats=AmpModul_Feats(1:end-1);                                         %% ---_Fix This 

                    Effort_e_Amp=(AmpModul_Feats-AmpModul_MAV_Noise)/(AmpModul_MAV_MAX-AmpModul_MAV_Noise)*100;
                    Effort_Amp=Effort_e_Amp+Effort_o;
                    Effort_Amp_Indv=Effort_e_Amp+Effort_o_Indv;

                    Effort_e_Amp_err=Effort_e_Amp-Effort_Fv_prime;
                    Effort_Amp_err=Effort_Amp-Effort_Fv_prime;
                    Effort_Amp_Indv_err=Effort_Amp_Indv-Effort_Fv_prime;

                    Force=S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames(:,KeepInd);

                    SNR_F= 20*log(RepTableMat(iTrial,:).('Target_Level')/std(Force(ConstInd)));
                    SE_F= std(Force(ConstInd))/sqrt(length(Force(ConstInd)));              

                    SNR_Amp= 20*log(mean(Effort_e_Amp(ConstInd))/std(Effort_e_Amp(ConstInd)));
                    SE_Amp= std(Effort_e_Amp)/sqrt(length(Effort_e_Amp));

                    ErrMat=[ErrMat; mean(Effort_e_MAV(ConstInd)) mean(Effort_MAV(ConstInd)) mean(Effort_MAV_Indv(ConstInd))...
                        mean(Effort_e_Amp(ConstInd)) mean(Effort_Amp(ConstInd)) mean(Effort_Amp_Indv(ConstInd))...
                        mean(Effort_e_MAV_err(ConstInd)) mean(Effort_MAV_err(ConstInd)) mean(Effort_MAV_Indv_err(ConstInd))...
                        mean(Effort_e_Amp_err(ConstInd)) mean(Effort_Amp_err(ConstInd)) mean(Effort_Amp_Indv_err(ConstInd))...
                        SNR_F SE_F SNR_Amp SE_Amp TypeLabel FiltLabel TrialLabel ExpLabel TestFolders(iTest) LogModel(iModel)...
                        EffortType S.(TestLabel).(ExpLabel).RepTableMat(iTrial,1:7) ];

%                     % 4- Saving the Results
    %                     S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst=table();
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_o')=array2table(Effort_o);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_o_Indv')=array2table(Effort_o_Indv);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_f')=array2table(Effort_f);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_Fv_prime')=array2table(Effort_Fv_prime);

                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_e_MAV')=array2table(Effort_e_MAV);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_MAV')=array2table(Effort_MAV);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_MAV_Indv')=array2table(Effort_MAV_Indv);

                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_e_Amp')=array2table(Effort_e_Amp);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_Amp')=array2table(Effort_Amp);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_Amp_Indv')=array2table(Effort_Amp_Indv);

                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_e_MAV_err')=array2table(Effort_e_MAV_err);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_MAV_err')=array2table(Effort_MAV_err);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_MAV_Indv_err')=array2table(Effort_MAV_Indv_err);

                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_e_Amp_err')=array2table(Effort_e_Amp_err);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_Amp_err')=array2table(Effort_Amp_err);
                    S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(:,'Effort_Amp_Indv_err')=array2table(Effort_Amp_Indv_err);
                end
            end
        end
    end
end

Effort_Est=array2table(ErrMat,'VariableNames',["Effort_e_MAV" "Effort_MAV" "Effort_MAV_Indv"...
    "Effort_e_Amp" "Effort_Amp" "Effort_Amp_Indv" "Effort_e_MAV_err" "Effort_MAV_err"...
    "Effort_MAV_Indv_err" "Effort_e_Amp_err" "Effort_Amp_err" "Effort_Amp_Indv_err" ...
    "SNR_Force" "SE_Force" "SNR_Amp" "SE_Amp" "Occ_Type" "Filt" "Trial" "Exp" "Test"...
    "LogModel" "EffortType" "Target_Level" "Stim_Force" "Voli_Force" "VoliMVC" "StimMVC" "PW" "Done"]);

writetable( Effort_Est,'occ_est_error4.csv')

%% Dropped Frames based Occlusion Estimation 
% Occ = 0-stim - avg(dropped frames)
% TestFolders=["jan7" "jan11" "apr20" "mar16"];

TimeRange=[10 15] ;
AnaLabel=sprintf("%s_ana",TestFolders{1});
ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');

% occ_table= readtable('occlusion_v5.csv');

FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);

% Occ=table([],[],[],[],[],[],[],[],[],'VariableNames',["MAV" "MAV_Type"...
%     "NormCoef" "Repeat" "Trial" "MVC_Stim" "MVC_Voli" "Test" "NumofDropped"]);

Occ=[];
NormMVC=40;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    StimRange=S.(TestLabel).(ExpLabel).StimConstantRange;
    % StimRange=TimeRange;
    FrameRange=[StimRange(1)*stim_freq StimRange(2)*stim_freq];
    FrameRangeInd=[FrameRange(1): FrameRange(2)];
    NumofDropped=S.(TestLabel).(ExpLabel).num_of_dropped;
    NumofTrial=S.(TestLabel).(ExpLabel).NumofTrials;

    StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;

    sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
    MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps');
    vMVCLevs=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');

    Ind=(sMVCzero==0);
    % Ind_Reps=(sMVCzero_Reps==0);
    p=polyfit(vMVCLevs(Ind),MAVMean(Ind),1);
    NormCoef=polyval(p,100); %% --- >>> Update based on mainanalysis
    % NormCoef=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;

    % NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX_theo
    NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX;

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    if MAV_MAXMethods(iTest)=="Fitted"
        NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;
        AmpModul_max_theo=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;
    else
        NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX;
        AmpModul_MAV_MAX=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;        
    end


    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    EffortLabel=S.(TestLabel).ExpPar.EffortType;    
    
    DroppedFiltLabel=DroppedFiltLabels(iTest);
    NoStimFiltLabel=NoStimFiltLabels(iTest);
   
    DropOccTest= [[] [] [] [] [] [] [] []];
    
    for iStim=2:length(StimMVCLevels)
        for iVoli=1:length(VoliMVCLevels)
            
            TrialNums=find_trialnum(VoliMVCLevels(iVoli),StimMVCLevels(iStim),... 
                S.(TestLabel).(ExpLabel).RepTableMat);

            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==NoStimFiltLabel;

            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('vMVC');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('sMVC');
            NoStim=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd & vMVC==VoliMVCLevels(iVoli) & sMVC==0,:).('MAV_Mean');
            
            for iRep=1:length(TrialNums)
                TrialLabel=sprintf("Trial_%d",TrialNums(iRep));
                
                DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;

                FiltInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat.('Filt')==DroppedFiltLabel;
                MAV_dropped_temp=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(FiltInd,:).('MAV');
                MAV_dropped= double(MAV_dropped_temp(FrameRange(1)<=DroppedFrames & DroppedFrames<=FrameRange(2),:));
                
                EffortMea=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortLabel);
                
                RowInd=OccTable.('Test')==TestFolders(iTest) & OccTable.('Trial')==string(TrialNums(iRep)) & ...
                    OccTable.('Feat')=="EffortMea" & OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau" & OccTable.('MVC_Stim')~="0";
               
                occEffort(iRep)=OccTable(RowInd,:).('Occ_mvc');
                occHybrid(iRep)=OccTable(RowInd,:).('Occ_Hybrid_mvc');
                occDropped2(iRep)=OccTable(RowInd,:).('Occ_Dropped');
                
                occDropped(iRep)=(mean(MAV_dropped)-mean(NoStim))/NormCoef*100;  %%----------> Use earlier calculations here as well
                NoStimMAV(iRep)=NoStim(iRep)/NormCoef*100;
                DroppedMAV(iRep)=mean(MAV_dropped)/NormCoef*100;
                Repeats(iRep)=sprintf("Rep_%d",iRep);
                Trials(iRep)=TrialLabel;
                StimMVC(iRep)=StimMVCLevels(iStim);
                VoliMVC(iRep)=VoliMVCLevels(iVoli);
                Tests(iRep)=string(TestFolders(iTest));
                Filters(iRep)=string(DroppedFiltLabel);
                NumDropped(iRep)=sprintf("%d_Drops",NumofDropped);
                NormCoefs(iRep)=NormCoef;
                
                DropOccTest=[DropOccTest; [occDropped(iRep) NoStimMAV(iRep) DroppedMAV(iRep) occEffort(iRep) occHybrid(iRep)]'...
                    ["OccDropped" "NoStim" "DroppedMAV" "OccEffort" "OccHybrid"]' ones(5,1)*NormCoefs(iRep) strings(5,1)+Repeats(iRep)...
                    strings(5,1)+Trials(iRep) ones(5,1)*StimMVC(iRep) ones(5,1)*VoliMVC(iRep) strings(5,1)+Filters(iRep) strings(5,1)+Tests(iRep)...
                    strings(5,1)+EffortLabel strings(5,1)+NumDropped(iRep) ];
                
                for iFilt=1:length(FiltLabels)
                    FiltLabel=FiltLabels(iFilt);

                    MAV_Filt=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameRange(1):FrameRange(2),:).('MAV_vEMG'))/NormCoef*100;

                    DropOccTest=[DropOccTest; MAV_Filt FiltLabel NormCoefs(iRep) Repeats(iRep) Trials(iRep) StimMVC(iRep) VoliMVC(iRep)...
                        DroppedFiltLabel Tests(iRep) EffortLabel NumDropped(iRep) ];
                end
            end
        end
    end
    
    DropOccTest=array2table(DropOccTest,'VariableNames',["MAV" "MAV_Type" "NormCoef" "Repeat" "Trial"...
        "MVC_Stim" "MVC_Voli" "Dropped_Filt" "Test" "EffortType" "NumofDropped"]);
    
    S.(AnaLabel).(ExpLabel).DropOccTest=DropOccTest;
    Occ=[Occ; DropOccTest];
end

writetable( Occ, 'dropped_occ6.csv')

%% Saving the results

% suffix="-func-";
% save_test(TestFolders,S,suffix)

%% Time Constants of RCCurve Trials

%To Do List
% 1- remove filtering of force^
% 2- indicing simplifications ^
% 3- Saving fixes: dont save the averages, save all the coefficients ^
% 4- Plotting presentations^
% 5- make this part of main analysis
% 6- EMG t
% # Filter Design

% d1 = designfilt("lowpassiir",'FilterOrder',3, ...
%     'HalfPowerFrequency',0.01,'DesignMethod',"butter"); %,'SampleRate',fs
% 
% % fvtool(d1)
% 
% LPPass = 30;
% LPStop = 100;
% Ap = .1;

% clc
% clear F_mat IndTau PWPoints ID_PW 
% Tau_stats=[];
% Tau=[];
% for iTest=1:length(TestFolders) 
%     TestLabel=sprintf("%s_test",TestFolders(iTest));
%     AnaLabel=sprintf("%s_ana",TestFolders(iTest));
% 
%     ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
%     ExpLabel=ExpLabels(str2double(S.(TestLabel).ExpPar.ExpTable(2,:).('RC')));
%     Tau_stats_test=[];
%     Tau_test=[];
%     if ~S.(TestLabel).(ExpLabel).ExpRun; continue; end
% 
%     NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
%     DataIndTable= S.(TestLabel).ExpPar.DataIndTable;
%     DataVars=string(DataIndTable.Properties.VariableNames);
%     iForce=(S.(TestLabel).ExpPar.DataIndTable.("Force"));
%     iTime=(S.(TestLabel).ExpPar.DataIndTable.("Time"));
% 
%     PWVal=S.(TestLabel).(ExpLabel).PWTrials;
%     RCVar=S.(TestLabel).(ExpLabel).RCVar;
%     MVC=S.(TestLabel).(S.(TestLabel).ExpPar.ExpTable(1,:).('MVC')).MVC;
%     TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
% 
%     AvgRangePostOff=[TurnOffTime+.5 TurnOffTime+.8];
%     AvgRangePreOff=[TurnOffTime-1.5 TurnOffTime]; 
% 
%     stim_freq=S.(TestLabel).ExpPar.stim_freq;
%     PostOffFrameRange=round([AvgRangePostOff(1)*stim_freq AvgRangePostOff(2)*stim_freq]);
%     PostOffFrameRangeInd=[PostOffFrameRange(1): PostOffFrameRange(2)];
%     PreOffFrameRange=round([AvgRangePreOff(1)*stim_freq AvgRangePreOff(2)*stim_freq]);
%     PreOffFrameRangeInd=[PreOffFrameRange(1): PreOffFrameRange(2)];
%     fs=S.(TestLabel).ExpPar.fs;
%     clear Taus
%     for iTrial=1:NumofTrials
%         TrialLabel=sprintf("Trial_%d",iTrial);
% 
%         F=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Force');
%         T=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
% 
%         F_filtered=F;
% 
%         Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
%         PreAvgForce= mean(F_filtered(Ind));
%         Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
%         PostAvgForce=mean(F_filtered(Ind));
% 
%         tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1))-PostAvgForce,-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce-PostAvgForce);
% 
%         IndTau=round((TurnOffTime+tau)*fs);
% 
%         MVCValofPWs=gompertz(RCVar,PWVal(iTrial))/MVC*100+0.00001;
% 
%         Tau_test=[Tau_test; tau IndTau PWVal(iTrial) MVCValofPWs PreAvgForce PostAvgForce iTrial TestFolders(iTest)  S.(AnaLabel).(ExpLabel).(TrialLabel).Redo ] ;
%         Tau=[Tau; Tau_test ] ;
% 
%         F_mat(:,1)=T';
%         VarNms(1)={'Time'};
%         F_mat(:,iTrial+1)=F_filtered()';
%         VarNms{iTrial+1}=char(TrialLabel);
% 
%     end
% 
% 
%     VarNames={'TimeConst','TauInd','PW','StimMVC','PreAvgForce','PostAvgForce','TrialNum','Test','Redo'};
%     Tau_table_test=array2table(Tau_test,'VariableNames',VarNames);
% 
% 
%     S.(AnaLabel).(ExpLabel).TimeCons.TurnOffTime=TurnOffTime;
%     S.(AnaLabel).(ExpLabel).TimeCons.AvgRangePostOff=AvgRangePostOff;
%     S.(AnaLabel).(ExpLabel).TimeCons.AvgRangePreOff=AvgRangePreOff;
%     S.(AnaLabel).(ExpLabel).TimeCons.Tau_test=Tau_test;
% 
%     F_table=array2table(F_mat,'VariableNames',VarNms);
%     clear F_mat
%     S.(AnaLabel).(ExpLabel).TimeCons.F_filtered=F_table;
%     S.(AnaLabel).(ExpLabel).Tau_table_test=Tau_table_test;
%     S.(AnaLabel).(ExpLabel).F_table=F_table;
% 
%     [ID_PW{iTest},PWPoints{iTest}]=findgroups(Tau_table_test(:,3));
%     for PWPointsInd=1:height(PWPoints{iTest})
% 
%         if isempty( PWPointsInd)
%            disp('PW value was not found') 
%            return
%         end
% 
%         PWInd(:,iTest)=ID_PW{iTest}== PWPointsInd;
%         TauVals=str2double(Tau_table_test(PWInd(:,iTest),:).('TimeConst'));
%         MVCVals=str2double(Tau_table_test(PWInd(:,iTest),:).('StimMVC'));
% 
%         Tau_stats_test=[Tau_stats_test; mean(TauVals)  std(TauVals) PWPoints{iTest}(PWPointsInd,:).('PW') mean(MVCVals) TestFolders(iTest)];
%         Tau_stats=[Tau_stats; Tau_stats_test];
%     end
%     VarNames2=["Mean" "Std" "PW" "StimMVC" "Test"];
%     Tau_stats_test=array2table(Tau_stats_test,'VariableNames',VarNames2);
% 
%     DirLabelCSV=sprintf('%s/%s_test_tau.csv',TestFolders(iTest),TestFolders(iTest));
%     writetable(Tau_table_test, DirLabelCSV)
%     DirLabelCSV=sprintf('%s/%s_taustats.csv',TestFolders(iTest),TestFolders(iTest));
%     writetable( Tau_stats_test, DirLabelCSV)
% end
% 
% Tau=array2table(Tau,'VariableNames',VarNames);
% writetable( Tau,'tau_estimates3.csv')
% 
% Tau_stats=array2table(Tau_stats,'VariableNames',VarNames2);
% writetable( Tau_stats,'tau_est_stats.csv')
% 
% 
% %% # Plotting 
% clc
% % PWPoints{1:length(TestFolders)}
% close all
% PlotPW = [130 120 120 100 100 105];
% for iTest=1:length(TestFolders)
%     TestLabel=sprintf('%s_test',TestFolders{iTest});
% 
%     Ind(iTest)=find(double(PWPoints{iTest}.('PW'))==PlotPW(iTest))
% 
%     if isempty( Ind{iTest})
%         sprintf('PW value was not found for %s',TestFolders{iTest})
%         return
%     end
% 
%     PWInd{iTest}=ID_PW{iTest} == Ind(iTest);
%     Tau{iTest}(PWInd{iTest},:);
%     TrialNums{iTest}=Tau{iTest}(:,'TrialNum');
% 
%     PWInd_F=logical([1 PWInd{iTest}']');
%     F_PW=K.(TestLabel).(ExpLabel).F_table(:,PWInd_F);
%     taus=Tau_table_test.('TimeConst')(PWInd{iTest});
%     trialnum=Tau_table_test.('TrialNum')(PWInd{iTest});
% 
%     %%
%     %Plotting here
%     h= figure(1);
%     set(h, 'Visible', 'on');
%     subplot(length(TestFolders),1,iTest)
%     plot(table2array(F_PW(:,1)),table2array(F_PW(:,2:end)),'LineWidth',2)
%     hold on
%     plot([0 5 6 10 10.00001 11],[0 0 PlotPW(iTest) PlotPW(iTest) 0  0]/20)
%     xlabel('Time(s)')
%     ylabel('Force(N)')
%     PWvalues=table2array(PWPoints{iTest});
%     TrialVals=table2array(TrialNums{iTest}(PWInd{iTest},'TrialNum'));
%     trials_str1=sprintf("Trials: " );
%     trials_str2=sprintf("%d, ",TrialVals );
%     PW_str=sprintf("PW: %d ",PlotPW(iTest) );
%     str_test=sprintf("Test: %s ",TestFolders{iTest} );
%     title(strcat(trials_str1,trials_str2,PW_str,",",str_test))
%     lgd_str{1}= sprintf ("Trial %d tau=%.2f",TrialVals(1),taus(1));
%     lgd_str{2}= sprintf ("Trial %d tau=%.2f",TrialVals(2),taus(2));
%     lgd_str{3}= sprintf ("Trial %d tau=%.2f",TrialVals(3),taus(3));
%     lgd_str{4}= sprintf ("Norm PW");
% 
%     legend(lgd_str,'Location','NorthEast','AutoUpdate','off')
%     xlim([5 11.9])
%     grid on
% 
% 
%     %% 
%     %annotations and mark important points
%     %%
% 
%     a = get(gca,'Children');
%     y1data = get(a, 'YData');
%     y1min=min( [min(y1data{1}) min(y1data{2}) min(y1data{3}) ]);
%     y1max=max( [max(y1data{1}) max(y1data{2}) max(y1data{3}) ]);
% 
%     PreAvgForce=table2array(Tau{iTest}(PWInd{iTest},'PreAvgForce'));
%     PostAvgForce=table2array(Tau{iTest}(PWInd{iTest},'PostAvgForce'));
%     AvgRangePreOff=K.(TestLabel).(ExpLabel).TimeCons.AvgRangePreOff;
%     AvgRangePostOff=K.(TestLabel).(ExpLabel).TimeCons.AvgRangePostOff;
%     TurnOffTime=K.(TestLabel).(ExpLabel).TimeCons.TurnOffTime;
% 
%     plot([TurnOffTime TurnOffTime]',[0 y1max]','--','Color','k')
%     plot([AvgRangePreOff(1) AvgRangePreOff(2)],[PreAvgForce PreAvgForce],'--','Color','k')
%     plot([AvgRangePostOff(1) AvgRangePostOff(2)],[PostAvgForce PostAvgForce],'--','Color','k')
% 
% end
end

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
xdft(abs(xdft)<=0.001)=0;
xdft = xdft(1:floor(length(x)/2)+1);
xdft(2:end-1) = 2*xdft(2:end-1);
psd = 1/(length(x)*fs)*abs(xdft).^2;
freq = 0:fs/length(x):fs/2;


MedFreq = medfreq(psd,freq);
MeanFreq = meanfreq(psd,freq);

end