%% mainana_func v2.0

%function S=ana_func(TestFolders,DroppedFiltLabels,NoStimFiltLabels,MAV_MAXMethods,TauTestsForce,TauTestsHand,AvgOcclusionTests,LogModel, BlankTime,gs_orders)


% %% Tidy Data
% 
TestFolders=["jan7" "jan11" "jan12" "aug22_24" "sep3_24" "sep4_24" "sep6_24"];
TestFolders=[ "sep4_24" "sep6_24" "oct17_24"];
TestFolders=["oct29_24", "oct31_24"];
TestFolders=["feb28_24", "feb29_24", "mar18_24", "mar20_24" ]; % these tests need manu_additoins
TestFolders=["oct29_24", "oct31_24", "nov14_24"];
TestFolders=["nov15_24" ];
TestFolders=[ "oct11" "oct18"];
TestFolders=["feb27"];
TestFolders=["oct25"];


for iTest=1:length(TestFolders)
    TestFolder=TestFolders(iTest);
    S=tidy_data(TestFolder);
    TestFile=sprintf("%s_test",TestFolder);
    str=sprintf('%s/%s',TestFolder,TestFile);
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
end


%% Defining Initial Parameters

% TestFolders=["mar20_24"];
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25" "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=[ "apr20" "oct11" "oct18"];
% TestFolders=["oct25"];
% TestFolders=[ "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=["jan7" "jan11" "jan12"  ];
% TestFolders=["jun20_24" "jul9_24" "jul21_24"  ]
% TestFolders=[ "jan7" "jan11" "jan12" "jun20_24" "jul9_24" "jul21_24" "jul31_24" ];
TestFolders=["feb28_24" "feb29_24" "mar18_24"  "mar20_24" ];
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24"];
TestFolders=[ "jan7" "jan11" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" "jan12" "aug22_24" "aug26_24"];
TestFolders=["aug29_24"];

TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24" "aug29_24"];
TestFolders=["jan7" "jan11" "jan12" "aug22_24" "sep3_24" "sep4_24" "sep6_24" "oct17_24" "oct18_24"];
TestFolders=["feb28_24" "feb29_24" "mar18_24"  "mar20_24" ];
% TestFolders=["jan7" "jan11" "jan12" "aug22_24" "sep3_24" "sep4_24" "sep6_24" "oct17_24" "oct18_24"];
% TestFolders=["sep3_24" "sep4_24" "oct17_24"];
TestFolders=["oct29_24", "oct31_24" "nov14_24"];
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24" "aug29_24" "sep3_24" "sep4_24" "sep6_24"];
% TestFolders=[ "jan7" "jan11" "jan12" "oct11" "nov15_24" "oct25"];
TestFolders=["oct29_24" "oct31_24" "nov14_24" ];
TestFolders=["sep6_24" "sep4_24" "aug26_24" "aug22_24" ];

NumofTests=length(TestFolders);
DroppedFiltLabels=strings(1,NumofTests)+"GS";
NoStimFiltLabels=strings(1,NumofTests)+"Unfilt";
BlankTime=ones(1,NumofTests)*0.004;
MAV_MAXMethods=strings(1,NumofTests)+"Real";
gs_orders=ones(1,NumofTests)*6;
ModelLog=false;

NumofTests=length(TestFolders);
DroppedFrameFilts=strings(1,NumofTests)+"GS";
NoStimFilts=strings(1,NumofTests)+"GS";
MAV_MAXMethods=strings(1,NumofTests)+"Real";

% AvgOcclusionTests=["sep6_24" "oct17_24" "oct18_24"];
AvgOcclusionTests=["jan7" "jan11" "jan12";
     "sep6_24" "sep4_24" "aug26_24" ];
% AvgOcclusionTests=["sep6_24" "sep4_24" "aug26_24"];

TauTestsForce=["jan7" "jan11" "jan12"]; % Average of 
TauTestsHand=["aug22_24" "nov15_24" "aug29_24"]; % Average of 

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
NumofTests=length(TestFolders);


% Data inject
textt="%s_test";
TestFiles=compose(textt,[TestFolders']);
S = load_test(TestFolders,TestFiles);

Tau_table= readtable('tau_estimates3.csv');
Tau_stats_table= readtable('tau_est_stats.csv');


%%Setting up initial parameters
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
%%Fixing PW invalid values issue (-inf)
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



%%Frames Matrix for Force and EMG
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

        if ExpLabel==ExpLabels(4)
            NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        end
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
            for iFrame=2:NumofFrames-1
                FrameLengths(iFrame)=-BegofFrames(iFrame)+BegofFrames(iFrame+1);
                FrameInd=BegofFrames(iFrame+1)+1-FrameLength:BegofFrames(iFrame+1);
                % # BP filt EMG
                x(:,iFrame)=BPFilt_EMG(FrameInd);
                % # Blanked EMG
                y(:,iFrame)=BlankEMG(FrameInd);
                % # Force
                f(:,iFrame)= Force(FrameInd);
                % # Time Frames
                t(:,iFrame)= Time(FrameInd);
                % # Trigger Frames
                tg(:,iFrame)= Trigger(FrameInd);
                    % # Hand
                h(:,iFrame)= Hand(FrameInd);
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

        NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
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

            DroppedFrames= ZeroInd(ZeroInd>=FrameRange(1) & ZeroInd<=FrameRange(2));
            DroppedEMG= S.(AnaLabel).(ExpLabel).(TrialLabel).BlankEMGFrames(:,DroppedFrames);

            S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames=DroppedFrames;
            S.(AnaLabel).(ExpLabel).(TrialLabel).ZeroInd=ZeroInd;
            S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG=DroppedEMG;
        end
    end
end


%%M-wave filtering  
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

%%EMG Features
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
                        % 
                        % [MedFreq_vEMG(iFrame), MeanFreq_vEMG(iFrame)]=MedMeanFreq(vEMGFrames(:,iFrame),fs);
                        % [MedFreq_MWaves(iFrame), MeanFreq_MWaves(iFrame)]=MedMeanFreq(MWaveFrames(:,iFrame),fs);
                        % 
                        % SSC_vEMG(iFrame)=NumSsc(vEMGFrames(:,iFrame));
                        % SSC_MWave(iFrame)=NumSsc(MWaveFrames(:,iFrame));
                        % 
                        % ZC_vEMG(iFrame)=NumZc(vEMGFrames(:,iFrame));
                        % ZC_MWave(iFrame)=NumZc(MWaveFrames(:,iFrame));

    %%%% Removed and replaced to improve run time of the code
                        % MedFreq_vEMG(iFrame)=medfreq(vEMGFrames(:,iFrame),fs);
                        % MedFreq_MWaves(iFrame)=medfreq(MWaveFrames(:,iFrame),fs);
                        % 
                        % MeanFreq_vEMG(iFrame)=meanfreq(vEMGFrames(:,iFrame),fs);
                        % MeanFreq_MWaves(iFrame)=meanfreq(MWaveFrames(:,iFrame),fs);
    %%%%    
                    MedFreq_vEMG(iFrame)=NaN;
                    MedFreq_MWaves(iFrame)=NaN;
                    MeanFreq_vEMG(iFrame)=NaN;
                    MeanFreq_MWaves(iFrame)=NaN;
    %%%%
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
                         table(MAV_vEMG',MAV_MWaves',MedFreq_vEMG',MeanFreq_vEMG',...
                         MedFreq_MWaves',MeanFreq_MWaves','VariableNames',...
                         ["MAV_vEMG" "MAV_MWave" "MedFreq_vEMG" "MeanFreq_vEMG"...
                         "MedFreq_MWave" "MeanFreq_MWave"]);
                    
                     if iExp ~= 1 % dont do this for MVC for having such small trial times less than the filter length
                         
                        FiltMAV_vEMG = filtfilt(lpFilt,MAV_vEMG);
                        FiltMedFreq_vEMG = filtfilt(lpFilt,MedFreq_vEMG);
                        FiltMeanFreq_vEMG = filtfilt(lpFilt,MeanFreq_vEMG);
                        FiltMAV_MWaves = filtfilt(lpFilt,MAV_MWaves);
                        FiltMedFreq_MWaves = filtfilt(lpFilt,MedFreq_MWaves);
                        FiltMeanFreq_MWaves = filtfilt(lpFilt,MeanFreq_MWaves);
                        % FiltSSC_vEMG = filtfilt(lpFilt,SSC_vEMG);
                        % FiltZC_vEMG = filtfilt(lpFilt,ZC_vEMG);
                        % FiltSSC_MWave = filtfilt(lpFilt,SSC_MWave);
                        % FiltZC_MWave = filtfilt(lpFilt,ZC_MWave);
    
                        Label=sprintf("FiltFeats%s",FieldSuffix(iField));
                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                            table(FiltMAV_vEMG',FiltMedFreq_vEMG',FiltMeanFreq_vEMG',...
                            FiltMAV_MWaves',FiltMedFreq_MWaves',FiltMeanFreq_MWaves',...
                            'VariableNames',[ "FiltMAV_vEMG" "FiltMedFreq_vEMG" "FiltMeanFreq_vEMG"...
                            "FiltMAV_MWave" "FiltMedFreq_MWave" "FiltMeanFreq_MWave"]);
    
                        FatMAV_vEMG = filtfilt(FatFilt,MAV_vEMG);
                        FatMedFreq_vEMG = filtfilt(FatFilt,MedFreq_vEMG);
                        FatMeanFreq_vEMG = filtfilt(FatFilt,MeanFreq_vEMG);
                        FatMAV_MWaves = filtfilt(FatFilt,MAV_MWaves);
                        FatMedFreq_MWaves = filtfilt(FatFilt,MedFreq_MWaves);
                        FatMeanFreq_MWaves = filtfilt(FatFilt,MeanFreq_MWaves);
                        % FatSSC_vEMG = filtfilt(FatFilt,SSC_vEMG);
                        % FatZC_vEMG = filtfilt(FatFilt,ZC_vEMG);
                        % FatSSC_MWave = filtfilt(FatFilt,SSC_MWave);
                        % FatZC_MWave = filtfilt(FatFilt,ZC_MWave);
    
                        Label=sprintf("FatFeats%s",FieldSuffix(iField));                    
                        S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(Label)=...
                            table(FatMAV_vEMG',FatMedFreq_vEMG',FatMeanFreq_vEMG',...
                            FatMAV_MWaves',FatMedFreq_MWaves',FatMeanFreq_MWaves',...
                            'VariableNames',[ "Filt_MAV_vEMG" "Filt_MedFreq_vEMG" "Filt_MeanFreq_vEMG"...
                            "Filt_MAV_MWave" "Filt_MedFreq_MWave" "Filt_MeanFreq_MWave" ]);
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

        TrialNum=S.(TestLabel).(ExpLabel).ListedNumofTrials;
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

RCMeanTime=[9 10]; % Calculating the means at time [8 10]
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
TrialStats=[];
TrialStatsAll=[];
Mean_TrialStatsAll=[];
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    VoliLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
    StimLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    TrialType=S.(TestLabel).(ExpLabel).TrialType;
    StimConstantRange=S.(TestLabel).(ExpLabel).StimConstantRange-2;
    MeanFrame=[StimConstantRange(1)*stim_freq StimConstantRange(2)*stim_freq];
    StimConstantRangeInd=[MeanFrame(1):MeanFrame(2)];

    VoliRange=StimConstantRange-2;
    VoliFrame=[VoliRange(1)*stim_freq VoliRange(2)*stim_freq];
    VoliRangeInd=[VoliFrame(1):VoliFrame(2)];
    
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
                        
                        % Mean_MAVDropped(TrialsInd(iTrial,c),1)=mean(DropFrames(MeanFrame(1)<=DroppedFrames & DroppedFrames<=MeanFrame(2))); 
                        % Std_MAVDropped(TrialsInd(iTrial,c),1)=std(DropFrames(MeanFrame(1)<=DroppedFrames & DroppedFrames<=MeanFrame(2)));    

                        DroppedFiltInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat.('Filt')==FiltLabel;
                        Mean_MAVDropped(TrialsInd(iTrial,c),1)=mean(str2double(S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(DroppedFiltInd,:).('MAV')));
                        Std_MAVDropped(TrialsInd(iTrial,c),1)=std(str2double(S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(DroppedFiltInd,:).('MAV')));
                    end
                    
                    MAV_Mean_reps(TrialsInd(iTrial,c),1)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(StimConstantRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),2)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(StimConstantRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),3)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(StimConstantRangeInd,:).('Amp_MAV_vEMG'));
                    MAV_Mean_reps(TrialsInd(iTrial,c),4)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(StimConstantRangeInd,:).('Amp_MAV_vEMG'));

                    MAV_Mean_Voli_reps(TrialsInd(iTrial,c),1)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(VoliRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_Voli_reps(TrialsInd(iTrial,c),2)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(VoliRangeInd,:).('MAV_vEMG'));
                    MAV_Mean_Voli_reps(TrialsInd(iTrial,c),3)=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(VoliRangeInd,:).('Amp_MAV_vEMG'));
                    MAV_Mean_Voli_reps(TrialsInd(iTrial,c),4)=std(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(VoliRangeInd,:).('Amp_MAV_vEMG')); 

                    MAV_Mean_reps(TrialsInd(iTrial,c),5)=Mean_MAVDropped(TrialsInd(iTrial,c),1);
                    MAV_Mean_reps(TrialsInd(iTrial,c),6)=Std_MAVDropped(TrialsInd(iTrial,c),1);
                    MAV_Mean_reps(TrialsInd(iTrial,c),7:9)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);
                    
                    Amp_Modul_Mean(TrialsInd(iTrial,c),3:5)=RepTableMat(TrialsInd(iTrial,c),[1,4,5]);
                    
                    TrialStats=[TrialStats; MAV_Mean_reps(TrialsInd(iTrial,c),1:4) MAV_Mean_Voli_reps(TrialsInd(iTrial,c),:) MAV_Mean_reps(TrialsInd(iTrial,c),5:9)...
                        FiltLabel TrialsInd(iTrial,c) TrialType ExpLabel TestLabel];                               
                end
                
                MAV_Mean(c,:) = [ c mean(MAV_Mean_reps(TrialsInd(:,c),1))  std(MAV_Mean_reps(TrialsInd(:,c),1)) mean(MAV_Mean_reps(TrialsInd(:,c),3))...
                    std(MAV_Mean_reps(TrialsInd(:,c),3)) mean(MAV_Mean_Voli_reps(TrialsInd(:,c),1))  std(MAV_Mean_Voli_reps(TrialsInd(:,c),1)) ...
                    mean(MAV_Mean_Voli_reps(TrialsInd(:,c),3)) std(MAV_Mean_Voli_reps(TrialsInd(:,c),3)) mean(Mean_MAVDropped(TrialsInd(:,c),1))...
                    std(Mean_MAVDropped(TrialsInd(:,c),1)) length(TrialsInd(:,c)) RepTableMat(TrialsInd(iTrial,c),[1 4 5 ])];
                
                Mean_TrialStats=[Mean_TrialStats; MAV_Mean(c,:) FiltLabel TrialsInd(iTrial,c) TrialType ExpLabel TestLabel];
            end
            
            c=c+1;
        end
    end
    
    % MAV_Mean_reps_table=table(double(TrialStats(:,1)),double(TrialStats(:,2)),double(TrialStats(:,3)),double(TrialStats(:,4)),double(TrialStats(:,5)),...
    %     double(TrialStats(:,6)),double(TrialStats(:,7)),double(TrialStats(:,8)),double(TrialStats(:,9)),TrialStats(:,10),TrialStats(:,11),TrialStats(:,12),...
    %     TrialStats(:,13),'VariableNames',["MAV_Mean" "MAV_Std" "Amp_Mean" "Amp_Std" "Dropped_Mean" "Dropped_Std" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "Exp" "Test"]);

    % MAV_Mean_table=table(double(Mean_TrialStats(:,1)),double(Mean_TrialStats(:,2)),double(Mean_TrialStats(:,3)),double(Mean_TrialStats(:,4)),double(Mean_TrialStats(:,5)),...
    %     double(Mean_TrialStats(:,6)),double(Mean_TrialStats(:,7)),double(Mean_TrialStats(:,8)),double(Mean_TrialStats(:,9)),double(Mean_TrialStats(:,10)),double(Mean_TrialStats(:,11)),...
    %     Mean_TrialStats(:,12),Mean_TrialStats(:,13),Mean_TrialStats(:,14),Mean_TrialStats(:,15),'VariableNames',["Unique_Trial" "MAV_Mean_Reps" "MAV_Std_Reps"...
    %     "Amp_Mean_Reps" "Amp_Std_Reps" "Dropped_Mean_Reps" "Dropped_Std_Reps" "NumofReps" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "Exp" "Test"]);

    MAV_Mean_reps_table=array2table(TrialStats,'VariableNames',["MAV_Mean" "MAV_Std" "Amp_Mean" "Amp_Std" "Voli_MAV_Mean" "Voli_MAV_Std" "Voli_Amp_Mean" "Voli_Amp_Std"...
        "Dropped_Mean" "Dropped_Std" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "TrialType" "Exp" "Test"]);
    
    MAV_Mean_table=array2table(Mean_TrialStats,'VariableNames',["Unique_Trial" "MAV_Mean_Reps" "MAV_Std_Reps" "Amp_Mean_Reps" "Amp_Std_Reps" "Voli_MAV_Mean_Reps"...
        "Voli_MAV_Std_Reps" "Voli_Amp_Mean_Reps" "Voli_Amp_Std_Reps" "Dropped_Mean_Reps" "Dropped_Std_Reps" "NumofReps" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial"...
        "TrialType" "Exp" "Test"]);
        
    S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table=MAV_Mean_reps_table;
    S.(AnaLabel).(ExpLabel).MAV_Mean_table=MAV_Mean_table;
    TrialStatsAll=[TrialStatsAll; TrialStats];
    Mean_TrialStatsAll=[Mean_TrialStatsAll; Mean_TrialStats];
    
end

TrialStatsAll_table=array2table(TrialStatsAll,'VariableNames',["MAV_Mean" "MAV_Std" "Amp_Mean" "Amp_Std" "Voli_MAV_Mean" "Voli_MAV_Std" "Voli_Amp_Mean" "Voli_Amp_Std"...
    "Dropped_Mean" "Dropped_Std" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "TrialType" "Exp" "Test"]);

Mean_TrialStatsAll_table=array2table(Mean_TrialStatsAll,'VariableNames',["Unique_Trial" "MAV_Mean_Reps" "MAV_Std_Reps" "Amp_Mean_Reps" "Amp_Std_Reps" "Voli_MAV_Mean_Reps"...
    "Voli_MAV_Std_Reps" "Voli_Amp_Mean_Reps" "Voli_Amp_Std_Reps" "Dropped_Mean_Reps" "Dropped_Std_Reps" "NumofReps" "TargetLevel" "vMVC" "sMVC" "Filt" "Trial" "TrialType"...
    "Exp" "Test"]);

% writetable(Mean_TrialStatsAll_table,"mean_trial_stats2.csv")
% writetable(TrialStatsAll_table,"trial_stats2.csv")

%%Normalize the Features over zero stim trials
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
            sMVC_Vec=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('sMVC'));
            vMVC_Vec=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('vMVC'));

            Voli_NormCoeff=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table(FiltInd,:).('MAV_Mean_Reps')...
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

%%Theoretical and Actual MVC MAV 

% NumofUpdate=4; % table mvc_table is updated 4 times inside loop below
% NumofVariables=5;
% mvc_table=strings(NumofUpdate*length(TestFolders),NumofVariables);
MVC_Percent=100;
sMVC_Ref=0;
AvgTimeMVC=2;
FiltLabel="Unfilt";
cm = lines(length(TestFolders));
mvc_table=[];
for iTest=1:length(TestFolders)
    
    % MAV_MAX = avg MAV at the last 2 seconds
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    EffortType=S.(TestLabel).ExpPar.EffortType;
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    stim_freq=S.(TestLabel).ExpPar.stim_freq;

    AvgInd=stim_freq*AvgTimeMVC-5;

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
    mvc_table=[mvc_table; [MAV_max "MAV" "Real" EffortType TestFolders(iTest) ]];
    mvc_table=[mvc_table; [AmpModul_MAV_max "Amp" "Real" EffortType TestFolders(iTest) ]];
    
    % mvc_table(iTest*NumofUpdate-3,:)=[MAV_max "MAV" "Real" EffortType TestFolders(iTest) ];        
    % mvc_table(iTest*NumofUpdate-2,:)=[AmpModul_MAV_max "Amp" "Real" EffortType TestFolders(iTest) ];
    
    % Theoretical MAV_MAX
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    sMVCLevs=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC'));
    vMVCLevs=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC'));
    
    MAVMean=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean'));
    Amp_Modul_Mean=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean'));
    
    Ind=(sMVCLevs==sMVC_Ref);
    
    poly1=polyfit(vMVCLevs(Ind),MAVMean(Ind),1);
    poly2=polyfit(vMVCLevs(Ind),Amp_Modul_Mean(Ind),1);
    
    MAV_max_theo=polyval(poly1,MVC_Percent);
    AmpModul_max_theo=polyval(poly2,MVC_Percent);
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    S.(AnaLabel).(ExpLabel).MAV_MAX_theo=MAV_max_theo;
    S.(AnaLabel).(ExpLabel).AmpModul_max_theo=AmpModul_max_theo;
    
    mvc_table=[mvc_table; [MAV_max_theo "MAV" "Theo" EffortType TestFolders(iTest) ]];
    mvc_table=[mvc_table; [AmpModul_max_theo "Amp" "Theo" EffortType TestFolders(iTest) ]];

    % mvc_table(iTest*NumofUpdate-1,:)=[ MAV_max_theo "MAV" "Theo" EffortType TestFolders(iTest) ];
    % mvc_table(iTest*NumofUpdate,:)=[ AmpModul_max_theo "Amp" "Theo" EffortType TestFolders(iTest) ];
    % 
    S.(AnaLabel).(ExpLabel).MVCTable= array2table(mvc_table,'VariableNames',["MVC_MAV" "MAV_Type" "Test_Type" "Effort_Type" "Test"]);


%%-----> Plotting the MVC_ MAV levels
    % 
    % figure(1) 
    % subplot(2,1,1)
    % plot(vMVCLevs(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    % hold on 
    % plot([unique(vMVCLevs(Ind)); 100],polyval(poly1,[unique(vMVCLevs(Ind)); 100]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest));
    % plot(100, MAV_max, '*','Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    % legend
    % xlabel('% Voli. MVC')
    % ylabel('Mean MAV')
    % % xlim([0 50])
    
end

% writetable(S.(AnaLabel).(ExpLabel).MVCTable,'mvc_table.csv')

%%Evaluate the filter accuracy
% NRMSE_s= sqrt(sum((MAV_voli-MAV_filt)^2)/T)/std(MAV_voli)
% NRMSE_v= sqrt(sum((MAV_voli_1-MAV_voli_2)^2)/T)/std(MAV_voli_1)
% TestFolders=["jan7" "jan11" "jan12" ];

c=1;

ExpstoAna= ["Occ"];
T=[];
% TestFolders=TestFolders(1);
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    TestName=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).("MVC");

    MVCTable=S.(AnaLabel).(ExpLabel).MVCTable;
    RowInd=MVCTable.("Test") == TestName & MVCTable.("MAV_Type") == "MAV" & MVCTable.("Test_Type") == "Real";
    MAV_max=str2double(MVCTable(RowInd,:).("MVC_MAV"))

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).("Occ");
    StimConstantRange=S.(TestLabel).(ExpLabel).StimConstantRange;
    StimRange=S.(TestLabel).(ExpLabel).StimRange;
    RampRange=[StimRange(1) StimConstantRange(1)];

    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameRange=[StimConstantRange(1)*stim_freq StimConstantRange(2)*stim_freq];
    FrameRangeInd=[FrameRange(1): FrameRange(2)];
    EffortType=S.(TestLabel).ExpPar.EffortType;

    NumofTrial=S.(TestLabel).(ExpLabel).ListedNumofTrials;
    RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat(:,1:6),...
        'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
        "MVC_Voli" "MVC_Stim" "PW"]);

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
            FiltLabel=FiltLabels(iFilt);

            MAV_filt=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameRangeInd,:).('MAV_vEMG');
            clear MAV_ref
            for iRefTrial=1:length(RefIndTrial)
                TrialLabel2=sprintf("Trial_%d",RefIndTrial(iRefTrial));
                
                MAV_ref=S.(AnaLabel).(ExpLabel).(TrialLabel2).Unfilt.Feats(FrameRangeInd,:).('MAV_vEMG');
                rmse(iRefTrial)=sqrt(sum((MAV_filt-MAV_ref).^2)/mean(MAV_ref))/length(MAV_filt);
                mav_diff(iRefTrial)=mean(MAV_filt-MAV_ref)/length(MAV_filt);
                Error=rmse(iRefTrial);
                
                Error_Type="RMSE";
                T=[T; Error Error_Type "false" MVC_Voli MVC_Stim iTrial RefIndTrial(iRefTrial) FiltLabel EffortType TestFolders(iTest)];

                MAV_ref_norm=MAV_ref/MAV_max;
                MAV_filt_norm=MAV_filt/MAV_max;
                rmse_norm(iRefTrial)=sqrt(sum((MAV_filt_norm-MAV_ref_norm).^2)/mean(MAV_ref_norm)) /length(MAV_filt);
                mav_diff_norm(iRefTrial)=mean(MAV_filt_norm-MAV_ref_norm)/length(MAV_filt);
                Error_norm=rmse_norm(iRefTrial);
                Error_Type="RMSE";
                T=[T; Error_norm Error_Type "true"  MVC_Voli MVC_Stim iTrial RefIndTrial(iRefTrial) FiltLabel EffortType TestFolders(iTest)];

                Error=mav_diff(iRefTrial);
                Error_Type="MAV_Diff";
                T=[T; Error Error_Type "false" MVC_Voli MVC_Stim iTrial RefIndTrial(iRefTrial) FiltLabel EffortType TestFolders(iTest)];
                
                Error=mav_diff_norm(iRefTrial);
                Error_Type="MAV_Diff";
                T=[T; Error Error_Type "true" MVC_Voli MVC_Stim iTrial RefIndTrial(iRefTrial) FiltLabel EffortType TestFolders(iTest)];                
            end
            
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).RMSE=mean(rmse);
            S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).mav_diff=mean(mav_diff);
            
        end  
    end
    % end
    RMSETable=array2table(T,'VariableNames',["Error",...
       "Error_Type" "Norm" "MVC_Voli" "MVC_Stim" "Trial" "Ref_Trial" "Filt" "Effort_Type" "Test"]);
    S.(AnaLabel).(ExpLabel).RMSETable=RMSETable;    
end

% writetable(RMSETable,'rmse3.csv')

%%
suffix="-jan3-";
save_test(TestFolders,S,suffix)

%% Occlusion Analysis
% # Determining Force, MAV Occlusion 
% Occ eqn F_o= Fs-Fs'   
% Fs : RC curve value with the corresponding PW -> R(PW)
% Fs': Force value at 3*tau

% To do :
% 1- incorporate redos: already incorporated earlier ^
% 2- a section goes to main analysis
% 3- what to do with mav target
% 4- Time const. for MAV drop
% 5- Combine other ways of calculating occ with this one

% AvgOcclusionTests=[  "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
clc
TauTests=["jan7" "jan11" "jan12"]; %% Tests for time constant calculation
% TauTests=["feb29_24" "mar18_24"  "mar20_24"]
TestInd=Tau_stats_table.('Test')==TauTests;

taus=Tau_stats_table(any(TestInd'),:).('Mean');
MeanTaus=mean(taus); % average time constant to be used 

OccRefs = [3*MeanTaus 4*MeanTaus 5*MeanTaus]; % referance time for occlusion to be calculated
TauLabels=["3*tau" "4*tau" "5*tau"];
OccType=["Occ_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc"];
NoiseAvgTime=[2 3];

Occ={};
for iTau=1:length(OccRefs(1))
    OccRef=OccRefs(iTau);

    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders(iTest));
        AnaLabel=sprintf("%s_ana",TestFolders(iTest));
        FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);

        EffortType=S.(TestLabel).ExpPar.EffortType;        
        DroppedFiltLabel=DroppedFrameFilts(iTest);
        NoStimFiltLabel=NoStimFilts(iTest);
        stim_freq=S.(TestLabel).ExpPar.stim_freq;        
        
        % MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;           %% MAV_MAX replacement 
        % AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;
       
        % MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX;
        % AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;
        
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
        if MAV_MAXMethods(iTest)=="Fitted"
            MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;           %% MAV_MAX replacement 
            AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;
        else
            MAV_max=S.(AnaLabel).(ExpLabel).MAV_MAX;
            AmpModul_MAV_max=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;
        end

        MVC=S.(TestLabel).(ExpLabel).MVC;
        
        % AvgTime=1.5;
        AvgInd=round(stim_freq*NoiseAvgTime(1):stim_freq*(NoiseAvgTime(2)));
        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        
        for iTrial=1:S.(TestLabel).(ExpLabel).ListedNumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial );

            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels(iFilt);
                MAV_Noise= mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats(AvgInd,:).('FiltMAV_vEMG'));
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MAV_Noise=MAV_Noise;

                AmpModul_MAV_Noise= mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(AvgInd,:).('Amp_MAV_vEMG'));
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModul_MAV_Noise=AmpModul_MAV_Noise;
            end
        end
        
        OccRefMargin= 0.2/stim_freq ; % in secs ------------> check if this is correct
        PreOffMargin=10/stim_freq; % Avg period before stim turning off 

        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
        NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat(:,1:7),...  %%-- Update this repmattable issue
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);
%         S.(TestLabel).(ExpLabel).RepTableMat=RepMatTable;       

        % FiltLabel="GS";
        % NoStimFiltLabel="GS";
        for iTrial=1:NumofTrials
            TrialLabel=sprintf("Trial_%d", iTrial);
            
            Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
            PostOffInd= TurnOffTime+OccRef-OccRefMargin<=Time & Time <=TurnOffTime+OccRef+OccRefMargin ;
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
            
            PostOffFrameInd=round((TurnOffTime+OccRef-OccRefMargin)*stim_freq:( TurnOffTime+OccRef+OccRefMargin)*stim_freq) ; %%%%%%%%%
            PreOffFrameInd=round((TurnOffTime-PreOffMargin)*stim_freq:( TurnOffTime)*stim_freq);
            % FiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('Filt')==FiltLabel;
            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('Filt')==NoStimFiltLabel;
            FiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==DroppedFiltLabel;
            TrialNums_reps=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Trial');

            MVC_Voli=RepMatTable(iTrial,:).('MVC_Voli');
            MVC_Stim=RepMatTable(iTrial,:).('MVC_Stim');
            sMVC=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'));
            vMVC=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC'));
            MAV_v=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table(NoStimFiltInd & MVC_Voli==vMVC & sMVC==0,:).('MAV_Mean_Reps')); % Are those means before or after the turn off ????
            AmpModul_MAV_v=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table(NoStimFiltInd & MVC_Voli==vMVC & sMVC==0,:).('Amp_Mean_Reps'));

            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==NoStimFiltLabel;
            DroppedMean_temp=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Dropped_Mean');

            DroppedMAV=double(DroppedMean_temp(iTrial==double(TrialNums_reps),:));
            vMVC_temp=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd,:).('vMVC');
            vMVC_iTrial=vMVC_temp(iTrial==double(TrialNums_reps),:);

            vMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('vMVC');
            sMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('sMVC');
            Mean_NoStimMAV=mean(double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd & vMVCTrials==vMVC_iTrial & double(sMVCTrials)==0,:).('MAV_Mean')));
            
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

                PreOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats_filtdropped(PreOffFrameInd,:).('MAV_vEMG'));
                PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats_filtdropped(PostOffFrameInd,:).('MAV_vEMG'));
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
                    .AmpModulFeats_filtdropped(PreOffFrameInd,:).('Amp_MAV_vEMG'));
                AmpModul_PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                    .AmpModulFeats_filtdropped(PostOffFrameInd,:).('Amp_MAV_vEMG'));

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

writetable( OccTable,'occlusion_v10.csv')
%% Occlusion Fitting  
%linear modeling for individual occ predictions 
    % Individual occlusion estimation
clc
clear lm_table
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20"];
lm_table=table();
CoefMat=[];
CoefMat2=[];
CoefMat3=[];

% AvgOcclusionTests=["feb28_24"]; % Generalized result will not match the 
% % indiv result for feb28_24 for this method because the slopes and averages will be slitly different
% AvgOcclusionTests=["jan7" "jan11" "jan12";
%                     "jun20_24" "jul9_24"	"jul21_24"]; 
% 
% AvgOcclusionTests=["jan7" "jan11" "jan12";
%      "aug22_24" "sep3_24" "sep4_24" ];
% AvgOcclusionTests=["sep6_24" "oct17_24" "oct18_24"];
clear  EffortTypeInd
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
        lm_table_half=table;        
        for iTest=1:length(EffortTestFolders)
            TestLabel=sprintf("%s_test",EffortTestFolders(iTest));
            lm_table=table();

            EffortType=S.(TestLabel).ExpPar.EffortType;
            RowInd=OccTable.('EffortType')==EffortType & OccTable.('Feat')=="EffortMea" & OccTable.('Filt')=="Unfilt"...
                & OccTable.('Tau')=="3*tau" & OccTable.('MVC_Stim')~="0" & OccTable.('Test')==EffortTestFolders(iTest);
            
% 1- Linear fitting for Individualized results
            lm_table.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
            lm_table.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
            lm_table.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
            lm_table.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
            lm_table.Test=categorical((OccTable(RowInd,:).('Test')));
            
            % lm_table.Test=reordercats(lm_table.Test, EffortTestFolders);
            mdl = fitlm(lm_table,'Occ~VoliMVC+StimMVC');
            CoefNum=table2array(mdl.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'))';
            CoefCat=table2array(mdl.Coefficients([1,4:end],'Estimate'));
            pVals=mdl.Coefficients.('pValue');

            % pValInd=[1 4:4+length(EffortTestFolders)-2]; % 1: p val for first test, 2,3: numerical variables, 4-:rest of the tests  
            % TestInd=lm_table.Test==EffortTestFolders(iTest);
            
            % if iTest==1
            %     Coefs(iTest,:)=[CoefCat(1) CoefNum];
            % else
            %     Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum];
            % end
            
            Coefs(iTest,:)=[CoefCat CoefNum];
            
            Effort_o=[];
            Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC+Coefs(iTest,3)*lm_table.StimMVC)];
            
            % CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) pVals(2) pVals(3) pVals(pValInd(iTest)) EffortTestFolders(iTest) EffortType...
            %     OccType(iType) boolean(1) boolean(0)], length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd)...
            %     lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]'];
             
            CoefMat3=[ CoefMat3; Coefs(iTest,:) pVals(1) pVals(2) pVals(3) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(0) ];

%%%%%%%% 2- Linear fitting for Half-Individualized results
            % RowInd=OccTable.('EffortType')==EffortType & OccTable.('Feat')=="EffortMea" & OccTable.('Filt')=="Unfilt"...
            %     & OccTable.('Tau')=="3*tau" & OccTable.('MVC_Stim')~="0";
            % 
            % lm_table_half.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
            % lm_table_half.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
            % lm_table_half.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
            % lm_table_half.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
            % lm_table_half.Test=categorical((OccTable(RowInd,:).('Test')));
            % 
            % lm_table_half.Test=reordercats(lm_table_half.Test,EffortTestFolders);
            % mdl = fitlm(lm_table_half,'Occ~VoliMVC+StimMVC+Test');
            % CoefNum=table2array(mdl.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'))';
            % CoefCat=table2array(mdl.Coefficients([1,4:end],'Estimate'));
            % pVals=mdl.Coefficients.('pValue');
            % pValInd=[1 4:4+length(EffortTestFolders)-2]; % 1: p val for first test, 2,3: numerical variables, 4-:rest of the tests  
            % TestInd=lm_table_half.Test==EffortTestFolders(iTest);
            % 
            % if iTest==1
            %     Coefs(iTest,:)=[CoefCat(1) CoefNum];
            % else
            %     Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum];
            % end
            % 
            % Effort_o=[];
            % Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table_half.VoliMVC(TestInd)+Coefs(iTest,3)*lm_table_half.StimMVC(TestInd))];
            % % 
            % % CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) pVals(2) pVals(3) pVals(pValInd(iTest)) EffortTestFolders(iTest) EffortType...
            % %     OccType(iType) boolean(1) boolean(0)], length(lm_table_half.Occ(TestInd)),1) Effort_o lm_table_half.Occ(TestInd)...
            % %     lm_table_half.LogOcc(TestInd) lm_table_half.VoliMVC(TestInd) lm_table_half.StimMVC(TestInd) [1:length(lm_table_half.Occ(TestInd))]'];
            % 
            % test_lbl=sprintf("%s_halfindiv",EffortTestFolders(iTest));
            % CoefMat3=[ CoefMat3; Coefs(iTest,:) pVals(pValInd(iTest)) pVals(2) pVals(3) test_lbl EffortType OccType(iType) boolean(1) boolean(0) ];
%%%%%%%%%%%%            
        end
        
        CoefMat=[ CoefMat; Coefs EffortTestFolders' strings(length(EffortTestFolders),1)+EffortType strings(length(EffortTestFolders),1)+OccType(iType) boolean(ones(length(EffortTestFolders),1))  boolean(zeros(length(EffortTestFolders),1))];
        
        mdlLog = fitlm(lm_table,'LogOcc~VoliMVC+StimMVC');
        CoefNum=table2array(mdlLog.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
        CoefCat=table2array(mdlLog.Coefficients([1,4:end],'Estimate'));
        pVals=mdlLog.Coefficients.('pValue');
        pValInd=[1 4:4+length(EffortTestFolders)-2]; % 1: p val for first test, 2,3: numerical variables, 4-:rest of the tests  
    %     Effort_o=[];
        for iTest=1:length(EffortTestFolders)
            % TestInd=lm_table.Test==EffortTestFolders(iTest);
    
            % if iTest==1
            %     Coefs(iTest,:)=[CoefCat(1) CoefNum'];
            % else
            %     Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum'];
            % end
            
            Effort_o=[];
            Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC+Coefs(iTest,3)*lm_table.StimMVC)];
            
            % CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) pVals(2) pVals(2) pVals(pValInd(iTest)) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(1)],...
            %     length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd)...
            %     lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];

            CoefMat3=[ CoefMat3; Coefs(iTest,:) pVals(1) pVals(2) pVals(3) EffortTestFolders(iTest) EffortType OccType(iType) boolean(1) boolean(1) ];
        end
        
        CoefMat=[ CoefMat; Coefs  EffortTestFolders' strings(length(EffortTestFolders),1)+EffortType strings(length(EffortTestFolders),1)+OccType(iType) boolean(ones(length(EffortTestFolders),1)) boolean(ones(length(EffortTestFolders),1))];
        

%----- 2- Linear fitting for Averaged results
        % lm_table.Test=string(lm_table.Test);
        RowInd=any(OccTable.('Test') == AvgOcclusionTests(iEffortType,:) & OccTable.('MVC_Stim')~="0",2);
        % lm_table_gen=lm_table(TestInd,:);
        lm_table_gen=table;
        lm_table_gen.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
        lm_table_gen.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
        lm_table_gen.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
        lm_table_gen.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
        lm_table_gen.Test=categorical((OccTable(RowInd,:).('Test')));
    
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
    
            % CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(0)],length(lm_table.Occ(TestInd)),1)...
            %     Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
            % 
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
    
            % CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(1)],length(lm_table.Occ(TestInd)),1)...
            %     Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
            
            CoefMat3=[ CoefMat3; CoefCat_gen CoefNum_gen' pVals' EffortTestFolders(iTest) EffortType OccType(iType) boolean(0) boolean(1)];
        end
    end
end

OccCoef=array2table(CoefMat,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "Test" "EffortType" "Type" "Indiv" "Log"]);

% OccCoef2=array2table(CoefMat2,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "pVal1" "pVal2" "pVal3" "Test" "EffortType" "Type" "Indiv" "Log"...
%     "Effort_o" "Occ" "LogOcc" "MVC_Voli" "MVC_Stim" "Trial"]);

OccCoef3=array2table(CoefMat3,'VariableNames',["Coeff1" "Coeff2" "Coeff3" "pVal1" "pVal2" "pVal3" "Test" "EffortType" "Type" "Indiv" "Log"]);
% writetable( OccCoef3,'occ_coef3.csv')

%% Effort Estimation
% Average Occlusion estimation
% OccCoef=[-1.7 -0.109 0.919];  % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fs_primeCoef=[2.117 0.061 0.118 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fv_primeCoef=[-1.361 0.891 0.848 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]

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

        % MAV_MAX=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;
        % AmpModul_max_theo=S.(AnaLabel).(ExpLabel).AmpModul_max_theo;

        MAV_MAX=S.(AnaLabel).(ExpLabel).MAV_MAX;
        AmpModul_MAV_MAX=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;        

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

        NumofTrials=S.(TestLabel).(ExpLabel).ListedNumofTrials;
        for iTrial=1:NumofTrials 
            TrialLabel=sprintf("Trial_%d",iTrial);

            DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
            PWofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames;
            StimMVC=S.(TestLabel).(ExpLabel).RepTableMat(iTrial,5);
            StimMVCofFrames=gompertz(RCVar,PWofFrames)/MVC*100+0.00001;

            KeepInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd(1:end-1);
            EffortTypeLabel=sprintf("%sFrames",EffortType);

            % if S.(TestLabel).ExpPar.EffortType=="Hand"
            %         Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).HandFrames(:,KeepInd))';
            %         Effort_f=Force;
            % 
            %     else
            %         Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames(:,KeepInd))';
            %         Effort_f=Force/F_MAX*100;
            % 
            % end

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
                
                % if EffortType=="Hand"
                %     Effort_o=Effort_o/100;
                %     Effort_o_Indv=Effort_o_Indv/100;
                % end

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

writetable( Effort_Est,'occ_est_error7.csv')

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
    FrameRangeInd=FrameRange(1): FrameRange(2);
    NumofDropped=S.(TestLabel).(ExpLabel).num_of_dropped;
    NumofTrial=S.(TestLabel).(ExpLabel).NumofTrials;

    StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;

    sMVCzero=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC')); 
    MAVMean=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps'));
    vMVCLevs=double(S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC'));
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');

    Ind=(sMVCzero==0);
    % Ind_Reps=(sMVCzero_Reps==0);
    p=polyfit(vMVCLevs(Ind),MAVMean(Ind),1);
    % NormCoef=polyval(p,100); %% --- >>> Update based on mainanalysis
    % NormCoef=S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX;

    % NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX_theo
    NormCoef=S.(AnaLabel).(ExpLabel).MAV_MAX;

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    EffortLabel=S.(TestLabel).ExpPar.EffortType;    
    
    % DroppedFiltLabel="GS";
    % NoStimFiltLabel="Unfilt";
    DroppedFiltLabel=DroppedFrameFilts(iTest);
    NoStimFiltLabel=NoStimFilts(iTest);
    DropOccTest= [[] [] [] [] [] [] [] []];
    
    for iStim=2:length(StimMVCLevels)
        for iVoli=1:length(VoliMVCLevels)
            
            TrialNums=find_trialnum(VoliMVCLevels(iVoli),StimMVCLevels(iStim),... 
                S.(TestLabel).(ExpLabel).RepTableMat);

            NoStimFiltInd=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt')==NoStimFiltLabel;

            vMVC=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('vMVC'));
            sMVC=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('sMVC'));
            NoStim=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(NoStimFiltInd & vMVC==VoliMVCLevels(iVoli) & sMVC==0,:).('MAV_Mean'));
            
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

% writetable( Occ, 'dropped_occ6.csv')

%%  Occlusion from dropped frames 
clc
MarginFromDropped=5;  % frames
TestStruct=sprintf("%s_test",TestFolders{iTest});

VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30];
stim_freq=S.(TestStruct).ExpPar.stim_freq;

boolean DroppedFrameInd;
TrialColor={'r', 'b', 'k'};
exp_lbl='Occ';
VarNames=["Frame_Val" "Norm_Val" "Target" "Frame" "Filt_Type" "MVC_Voli"...
    "MVC_Stim"  "Test"  "Feat" "Repeat" "Trial" "Drp_Order"];

NormVoliMVC=20;
g=strings(0,length(VarNames)-4);
x=ones(0,4);

g_dropped=strings(0,length(VarNames)-4);
x_dropped=ones(0,4);
x_feat=[];

AnaStruct=sprintf("%s_ana",TestFolders(1));

FeatLabels=string(S.(AnaStruct).AnaPar.FeatLabels);
ExpLabel=S.(AnaStruct).AnaPar.ExpTable(1,:).(exp_lbl);
sMVCNormLevel=0; % %MVC value for normalization coefficient
clear DrpOrderCat
for iFeat=1:1
    FeatLabel=FeatLabels(iFeat);
    
    for iTest=1:length(TestFolders)
        TestStruct=sprintf("%s_test",TestFolders{iTest});
        AnaStruct=sprintf("%s_ana",TestFolders{iTest});

        RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
        StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
        VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;

        for iVoli=1:length(VoliMVCLevels)
            sMVC_Vec=S.(AnaStruct).(ExpLabel).MAV_Mean_table.('sMVC');
            vMVC_Vec=S.(AnaStruct).(ExpLabel).MAV_Mean_table.('vMVC');

            Voli_NormCoeff=S.(AnaStruct).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps')...
                (sMVC_Vec==10 & vMVC_Vec==VoliMVCLevels(iVoli));

            for iStim=1:length(StimMVCLevels)
                sMVC=StimMVCLevels(iStim);
                vMVC=VoliMVCLevels(iVoli);
                IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

                for iRep=1:length(IndTrials)
                    TrialLabel=sprintf('Trial_%d',IndTrials(iRep));
                    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;

                    for iFilt=1:length(S.(AnaStruct).AnaPar.FiltLabels)
                        FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};
                        vEMGLabel=sprintf('%s_vEMG',FeatLabel);
                        NormLabel=sprintf('Norm_%s',vEMGLabel);
                        
                        DroppedFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(FeatLabel);
                        DroppedFrameNum=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;                        
                        
                        Feat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel);
                        Target=mean(S.(AnaStruct).(ExpLabel).TargetFrames)';
                        FrameNum=[1:length(Feat)];

                        x_feat=Feat;
                        x_target= Target(setdiff([1:length(Feat)+length(DroppedFeat)]',DroppedFrameNum));   
                        x_norm=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(NormLabel);
                        x_framenum=FrameNum';

                        lg=length(x_feat);
                        g_filt=strings(1,lg);
                        g_voli=strings(1,lg);
                        g_stim=strings(1,lg);
                        g_test=strings(1,lg);
                        g_feattype=strings(1,lg);
                        g_rep=strings(1,lg);
                        g_trial=strings(1,lg);
                        g_order=strings(1,lg);

                        g_filt(:)=string(FiltLabel);
                        g_voli(:)=sprintf("Voli_%d%%",(VoliMVCLevels(iVoli)));
                        g_stim(:)=sprintf("Stim_%d%%",(StimMVCLevels(iStim)));
                        g_test(:)=string(TestFolders(iTest));
                        g_feattype(:)=string(FeatLabel);
                        g_rep(:)=sprintf("Rep_%d",iRep);
                        g_trial(:)=string(TrialLabel);
                        g_order(:)="NonDropped";
                        
                        g=[g; g_filt' g_voli' g_stim' g_test' g_feattype' g_rep' g_trial' g_order' ];
                        x=[x ;x_feat x_norm x_target x_framenum ];
                        
                    end
                    
                    num_of_dropped=S.(TestStruct).(ExpLabel).num_of_dropped;
                    DroppedFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(FeatLabel);
                    DroppedFrameNum=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                    DroppedOrder=repmat([1:num_of_dropped],1,length(DroppedFeat)/num_of_dropped);
                    
                    x_droppedfeat=DroppedFeat;
                    x_droppedtarget=Target(DroppedFrameNum);
                    x_droppednorm=x_droppedfeat/Voli_NormCoeff; % ----------> this should go to main_analysis
                    x_droppedframenum=DroppedFrameNum;
                    
                    lg_dropped=length(x_droppedfeat);
                    
                    g_droppedfilt=strings(1,lg_dropped);
                    g_droppedvoli=strings(1,lg_dropped);
                    g_droppedstim=strings(1,lg_dropped);
                    g_droppedtest=strings(1,lg_dropped);                    
                    g_droppedfeattype=strings(1,lg_dropped);
                    g_droppedrep=strings(1,lg_dropped);
                    g_droppedtrial=strings(1,lg_dropped);
                    g_droppedorder=strings(1,lg_dropped);

                    g_droppedfilt(:)="Dropped";
                    g_droppedvoli(:)=sprintf("Voli_%d%%",(VoliMVCLevels(iVoli)));
                    g_droppedstim(:)=sprintf("Stim_%d%%",(StimMVCLevels(iStim)));
                    g_droppedtest(:)=string(TestFolders(iTest));
                    g_droppedfeattype(:)=string(FeatLabel);
                    g_droppedrep(:)=sprintf("Rep_%d",iRep);
                    g_droppedtrial(:)=string(TrialLabel);
                    g_droppedorder(:)=string(DroppedOrder);

                    g_dropped=[g_dropped; g_droppedfilt' g_droppedvoli' g_droppedstim' g_droppedtest' g_droppedfeattype' g_droppedrep' g_droppedtrial' g_droppedorder' ];
                    x_dropped=[x_dropped; x_droppedfeat x_droppednorm x_droppedtarget x_droppedframenum ];
                end
            end
        end
    end
end

Dropped_stats=array2table([[x ;x_dropped] [g; g_dropped] ],'VariableNames',[VarNames]);

% S.(AnaStruct).(ExpLabel).Dropped_stats=Dropped_stats;
% writetable( Dropped_stats, 'dropped_stats3.csv')


%% Plots before analysis
clc
close all
Trials=[ 1: 10]; 
TimeRange=[1 16];  % in seconds
% TestLabel=sprintf('%s_test',FolderNames);
ylims=5;
ExpLbl="Occ";
for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    TestName=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(ExpLbl);

    for iTrial=1:length(Trials)
        TrialLabel=sprintf('Trial_%d',Trials(iTrial));
        DataIndTable= S.(TestLabel).ExpPar.DataIndTable;

        Time=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,DataIndTable.("Time"));
        TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);
        EffortType=S.(TestLabel).ExpPar.EffortType;
        EMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("EMG"));
        Trigger=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("Trigger"));
        Force=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.(EffortType));
        PW=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,DataIndTable.("PW"));

        f=figure(iTest);
        f.Position = [100 100 1400 800];        
        subplot(2,length(Trials),iTrial)
        plot(Time(TimeInd),Trigger/20,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),EMG,'b','LineWidth',2)
        legend({'Trigger(a.u.)','Raw EMG (mV)'})
        title(sprintf('Test:%s, Trial: %d', TestName, Trials(iTrial)))
        xlabel('Time (s)')
        ylim(0.4*[-1 1])

        subplot(2,length(Trials),length(Trials)+iTrial)
        plot(Time(TimeInd),PW/200,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),Force,'b','LineWidth',2)
        legend({'PW (a.u.)','Meas. Effort'})
        title('Pulse Width and Force Signal')
        xlabel('Time (s)')
        ylim([-100 800])
    end

RepTableMat=S.(TestLabel).OccTrials.RepTableMat;
RepTableMat=[[1:length(RepTableMat)]' RepTableMat]

end

%% Plots after filtering and preproccess 
clc
close all
Trials=[ 1: 10]; 
TimeRange=[1 12];

% TestLabel=sprintf('%s_test',FolderNames);
ylims=5;
ExpNum=4;
for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    TestName=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpPar.ExpLabels(ExpNum);
    EffortType=S.(TestLabel).ExpPar.EffortType;

    for iTrial=1:length(Trials)
        TrialLabel=sprintf('Trial_%d',Trials(iTrial));
        DataIndTable= S.(TestLabel).ExpPar.DataIndTable;
        % TimeRange=S.(TestLabel).(ExpLabel).StimRange;
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
        TimeInd= Time>=TimeRange(1) & Time<=TimeRange(2);

        % RawEMG=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('EMG');
        BlankedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('BlankedEMG');

        Trigger=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('Trigger');
        EffortMeas=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).(EffortType);
        PW=S.(AnaLabel).(ExpLabel).(TrialLabel).data(TimeInd,:).('PW');

        f=figure(iTest);
        f.Position = [100 100 1400 800];                
        subplot(2,length(Trials),iTrial)
        plot(Time(TimeInd),Trigger/2000,'r','LineWidth',1)
        hold
        plot(Time(TimeInd),BlankedEMG,'b','LineWidth',2)
        legend({'Trigger(a.u.)','Raw EMG (mV)'})
        title(sprintf('Test:%s, Trial: %d', TestName, Trials(iTrial)))
        xlabel('Time (s)')
        % ylim([-1 1]*10^-3)
    end
OccTable=S.(TestLabel).OccTrials.RepTableMat;
OccTable=[[1:length(OccTable)]' OccTable]

end

%% Plotting, Debugging after Filtering
close all
Lbl='Occ';
PlotTrial=[ 10 11];
PlotTime=[ 12 13];
% PlotFrame=floor([ stim_freq*PlotTime(1) stim_freq*PlotTime(2)]);
PlotFrame=floor([ 278 282]);

TeststoPlot=TestFolders(end);

for iTest=1:length(TeststoPlot)
    TestLabel=sprintf("%s_test",TeststoPlot(iTest));
    AnaLabel=sprintf("%s_ana",TeststoPlot(iTest));
    DataLabels=S.(AnaLabel).AnaPar.DataLabels;
    TestName=S.(AnaLabel).AnaPar.TestName;
    ExpLabel=string(S.(TestLabel).ExpPar.ExpTable(1,:).(Lbl));

    % PlotFrame=floor([ 300 305]); % first frame is skipped
    % ExpLabels=S.(TestLabel).ExpPar.ExpLabels;
    % ExpstoAna=ExpLabels([2,3,4,5]);
    FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;

    clear EMGFrames MWaveFrames vEMG MWave MWaveFrames_filtdropped vEMGFrames_filtdropped vEMG_filtdropped MWave_filtdropped
    % vEMGFrames_filtdropped=zeros(FrameLength,PlotFrame(2)-PlotFrame(1)+1,length(FiltLabels));
    for iTrial=PlotTrial(1):PlotTrial(2)
        TrialLabel=sprintf('Trial_%d',iTrial);
        FrameLength=S.(AnaLabel).(ExpLabel).(TrialLabel).FrameLength;

        for iFilt=1:length(FiltLabels)
            FiltLabel=FiltLabels{iFilt};

            EMGFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));

            vEMG(:,iFilt)=reshape(EMGFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave(:,iFilt)=reshape(MWaveFrames(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

            vEMGFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            MWaveFrames_filtdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_filtdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));

            vEMG_filtdropped(:,iFilt)=reshape(vEMGFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            MWave_filtdropped(:,iFilt)=reshape(MWaveFrames_filtdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

            % vEMGFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            % MWaveFrames_withdropped(:,:,iFilt)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MWaveFrames_withdropped(1:FrameLength,PlotFrame(1):PlotFrame(2));
            % 
            % vEMG_withdropped(:,iFilt)=reshape(vEMGFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
            % MWave_withdropped(:,iFilt)=reshape(MWaveFrames_withdropped(:,:,iFilt),1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        end

        TimeFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TmFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));
        TriggerFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).TriggerFrames(1:FrameLength,PlotFrame(1):PlotFrame(2));

        Time=reshape(TimeFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));
        Trig=reshape(TriggerFrames,1,FrameLength*length(PlotFrame(1):PlotFrame(2)));

        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{2},TeststoPlot(iTest)),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,Trig/10000,'b','LineWidth',2)
        hold on
        plot(Time,vEMG(:,1),'k','LineWidth',2)

        plot(Time,vEMG(:,2),'r','LineWidth',2)
        legend({'Trigger(a.u)', 'Unfiltered', 'Comb vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d'...
            ,DataLabels{8},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,2),'r','LineWidth',2)
        legend({'Unfiltered', 'Comb MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{9},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

        figure('Name',sprintf('%s Results, Test: %s',FiltLabels{3},TeststoPlot(iTest)),'NumberTitle','off')
        subplot(2,1,1)
        plot(Time,Trig/10000,'b','LineWidth',2)
        hold on
        plot(Time,vEMG(:,1),'k','LineWidth',2)

        plot(Time,vEMG(:,3),'r','LineWidth',2)
        legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
        ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

        subplot(2,1,2)
        plot(Time,vEMG(:,1),'k','LineWidth',2)
        hold
        plot(Time,MWave(:,3),'r','LineWidth',2)
        legend({'Unfiltered', 'GS MWaves'})
        ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
            DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        title(ttl);
        xlabel(DataLabels{5})
        ylabel(DataLabels{1})

        % figure(10)
        % subplot(2,1,1)
        % plot(Time,Trig/10000,'b','LineWidth',2)
        % hold on
        % plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        % 
        % plot(Time,vEMG_filtdropped(:,3),'r','LineWidth',2)
        % legend({'Trigger(a.u)', 'Unfiltered', 'GS vEMG'})
        % ttl=sprintf('%s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
        %     DataLabels{10},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        % title(ttl);
        % xlabel(DataLabels{5})
        % ylabel(DataLabels{1})
        % 
        % subplot(2,1,2)
        % plot(Time,vEMG_filtdropped(:,1),'k','LineWidth',2)
        % hold
        % plot(Time,MWave_filtdropped(:,3),'r','LineWidth',2)
        % legend({'Unfiltered', 'GS MWaves'})
        % ttl=sprintf(' %s at %s, Test: %s, TrialNum: %d, Frame: %d-%d',...
        %     DataLabels{11},ExpLabel,TestName,iTrial,PlotFrame(1),PlotFrame(2)+1);
        % title(ttl);
        % xlabel(DataLabels{5})
        % ylabel(DataLabels{1})

    end
end

%% Plotting the dropped frames 
close all
iTrial=12;
TrialLabel=sprintf('Trial_%d',iTrial);
ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
FiltLabel="GS";
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    if ~S.(TestStruct).(ExpLabel).dropped 
                    figure(iTrial)
            subplot(length(TestFolders),1,iTest)
            text(0.3,0.5,'No Dropped Frames for This Trial')
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
        continue;
    else
       DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;   
        if isempty(DroppedFrames)
            figure(iTrial)
            subplot(length(TestFolders),1,iTest)
            text(0.3,0.5,'No Dropped Frames for This Trial')
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            continue;
        end

        DataLabels=S.(AnaLabel).AnaPar.DataLabels;
        DroppedEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedEMG(1:end-3,:);
        FiltDropped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).vEMGFrames_filtdropped(:,DroppedFrames);
        clear lgd

        for iFrame=1:length(DroppedFrames)
            figure(iTrial)
            subplot(length(TestFolders),2,iTest*2-1)
            plot(DroppedEMG(:,iFrame),'LineWidth',2)
            hold on
            lgd(iFrame)=sprintf("Frame %d",DroppedFrames(iFrame));
            legend(lgd,'Location','NorthWest');
            ttl=sprintf('Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            xlabel('Samples')
            ylabel('BP Filtered EMG')
            subplot(length(TestFolders),2,iTest*2)
            plot(FiltDropped(:,iFrame),'LineWidth',2)
            hold on 
            ttl=sprintf('Filtered Dropped Frames at Trial %d of %s, Test: %s',iTrial,ExpLabel,TestFolders{iTest});
            title(ttl);
            ylabel('BP Filtered and m-Wave Filtered')


            % ylim([3 -3]*10^-3)

        end
    end
end


%% Plotting EMG Features 
clc
AnaLabel=sprintf('%s_ana',TestFolders(1));
ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
TimeRange=[1 10];
TrialNum=[ 1: 12];
cm=lines(3);
FiltLabel="GS";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=ceil([stim_freq*TimeRange(1): stim_freq*TimeRange(2)]);

    for iTrial=1:length(TrialNum)

        TrialLabel=sprintf('Trial_%d',TrialNum(iTrial));
        AmpMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Amp_MAV_vEMG');
        AmpClipped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Clip_MAV_vEMG');   

        figure(iTrial)
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,AmpMAV,'LineWidth',2,'DisplayName',sprintf('AmpMAV, Trial: %d',TrialNum(iTrial)),'Color',cm(1,:))
        hold on 
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,AmpClipped,'DisplayName',sprintf('AmpClipped, Trial: %d',TrialNum(iTrial)),'Color',cm(2,:))
        title(TestFolders(iTest))
    end
%     title(ttl)
    legend('Location','NorthWest')
    grid on
end


%% Mean MAV Plotting
cm = lines(length(TestFolders));
sMVC=0;
FiltLabel="Unfilt";
Extrapolate=100;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    FiltInd=FiltLabel==S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');

    sMVCzero_Reps=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('sMVC'));
    MAVMean_Reps=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('MAV_Mean'));
    vMVCVal_Reps=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('vMVC'));

    % sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
    % MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps');
    % vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
    % Ind=(sMVCzero==sMVC);

    Amp_Modul_Mean=double(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(FiltInd,:).('Amp_Mean'));

    Ind_Reps=(sMVCzero_Reps==sMVC);
    p1=polyfit(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),1);

    p2=polyfit(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),1);

    figure(1)
    subplot(2,1,1)
    plot(vMVCVal_Reps(Ind_Reps),MAVMean_Reps(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p1,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    xlabel('% Voli. MVC')
    ylabel('Mean MAV')
    % xlim([0 50])

    subplot(2,1,2)
    plot(vMVCVal_Reps(Ind_Reps),Amp_Modul_Mean(Ind_Reps),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    hold on
    plot([vMVCVal_Reps(Ind_Reps); Extrapolate],polyval(p2,[vMVCVal_Reps(Ind_Reps); Extrapolate]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    ylabel('% Amp modul')
    xlabel('% Voli. MVC')
    % xlim([0 50])
end


%% Comparing Theoretical and Actual MVC MAV 

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
    
    AmpModul_MAV_max=max(TrialsAmp_MAV);

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
%%-----> Plotting the MVC_ MAV levels

    figure(1) 
    subplot(2,1,1)
    plot(vMVCLevs(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    hold on 
    plot([unique(vMVCLevs(Ind)); 100],polyval(poly1,[unique(vMVCLevs(Ind)); 100]),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest));
    plot(100, MAV_max, '*','Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    legend
    xlabel('% Voli. MVC')
    ylabel('Mean MAV')
    % xlim([0 50])
    
end


%% Plotting Adaptive Filter 
PlotRange=[ 5 10];
vMVC=10;
sMVC=0;
FiltLabel="GS";
% TestFolders="apr20";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    
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
    end
end



%% Comparing MWave and Voli MAV
% TestFolders= "feb29_24";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC');
    OccMeans=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(sMVC==0,:).('MAV_Mean');    
    OccStds=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(sMVC==0,:).('MAV_Std');
    MVCRange=S.(AnaLabel).(ExpLabel).MAV_Mean_table(sMVC==0,:).('vMVC');
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('RC');
    RCMeans=double(S.(AnaLabel).(ExpLabel).MAV_Mean.('Mean'));
    RCStds=double(S.(AnaLabel).(ExpLabel).MAV_Mean.('Std'));
    PWPoints=double(S.(AnaLabel).(ExpLabel).MAV_Mean.('PW'));
    PWPoints(end)=PWPoints(end);
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

%% Plotting NRMSE
clc
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders(iTest));
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    
    RMSETable=S.(AnaLabel).(ExpLabel).RMSETable;

    VoliLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
    StimLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    for iVoli=1:length(VoliLevels)
        VoliLevel=VoliLevels(iVoli);
    
        for iStim=1:length(StimLevels)
            StimLevel=StimLevels(iStim);
            
            NRMSE=zeros(0,length(FiltLabels));
            RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
            % stim_freq=S.(TestLabel).ExpPar.stim_freq;
            TrialsInd= find_trialnum(VoliLevel, StimLevel, RepTableMat);

            % RowInd=RMSETable
            for iTrial=1:length(TrialsInd)
                TrialLabel=sprintf("Trial_%d",TrialsInd(iTrial));
                
                NRMSE=[NRMSE; 
                    [S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabels(1)).RMSE'
                     S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabels(2)).RMSE'
                     S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabels(3)).RMSE']'  ];
            end
    
            filt=3;
            if StimLevel==0
                FiltLabelsPlot=FiltLabels(filt);
                % NRMSE=NRMSE(NRMSE~=0);
                NRMSE(:,1)=NRMSE(:,filt);
    
                NRMSE(:,3)=[];
                NRMSE(:,2)=[];            
            else
                FiltLabelsPlot=FiltLabels;
            end
            
            figure(2)
            subplot(length(StimLevels),length(VoliLevels),(iVoli-1)*(length(VoliLevels))+iStim)
            bar(FiltLabelsPlot,mean(NRMSE,1))
            [r,c]=size(NRMSE);
            hold 
            er = errorbar(categorical(FiltLabelsPlot),mean(NRMSE,1),-std(NRMSE,1),+std(NRMSE,1));    
            % er.Color = [0, 0, 0];
            % er.LineStyle = 'none'; 
            ylim([ 0 0.3])
            grid
            subplot(length(StimLevels),length(VoliLevels),iStim)
            title(sprintf("Stim: %d%%, n= %d(each bar)",StimLevel,r),'fontweight','bold','fontsize',14)
    
        end
        subplot(length(StimLevels),length(VoliLevels),(iVoli-1)*length(StimLevels)+1)
        ylabel(sprintf('Voli: %d%%',VoliLevel),'fontweight','bold','fontsize',14)
    
    end
end
%% Occ Trials Superimposed

close all
clc
lbl='Occ';
PlotVoli=3;
PlotStim=2;
cm=lines(length(TestFolders));
TimeRange=[1 17];
FiltLabel="GS";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    TestName=TestFolders(iTest);
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    MVC=S.(TestLabel).(ExpLabel).MVC;

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(lbl);
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    StimMVCLevels=S.(TestLabel).(ExpLabel).StimMVCVec;
    VoliMVCLevels=S.(TestLabel).(ExpLabel).VoliMVCVec;
    sMVC=StimMVCLevels(PlotStim);
    vMVC=VoliMVCLevels(PlotVoli);
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

    StimOff=S.(TestLabel).(ExpLabel).StimProfile;
    PlotRange=[ StimOff-10 StimOff+1];
    PreStimOffRange= [StimOff-1 StimOff ];
    PostStimOffRange=[StimOff+.2 StimOff+.3];
    ttl=sprintf('Stim: %d%%, Voli: %d%%',sMVC,vMVC);
    EffortLabel=S.(TestLabel).ExpPar.EffortType;
    for iTrial=1:length(IndTrials)
    IndTrials(iTrial)
        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        EffortMea=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortLabel);
        MAV=S.(AnaLabel).(ExpLabel).(TrialLabel).Unfilt.Feats.('MAV_vEMG');

        FiltMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FiltFeats.('FiltMAV_vEMG');
        AmpModulMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Amp_MAV_vEMG');

        MAVInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd;
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time'); 
        TimeInd=Time>TimeRange(1) & Time<TimeRange(2);

        DroppedMAV=[];
        DroppedInd=[];
        DroppedFeatFilt=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat.('Filt');
        if sMVC~=0
            DroppedMAV=str2double(S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(DroppedFeatFilt==FiltLabel,:).('MAV'));
            DroppedInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
        end
        f=figure(1);
        f.Position = [100 300 2000 800];    
        subplot(2,1,1)
        plot(Time(TimeInd),EffortMea(TimeInd),'Color',cm(iTest,:),'DisplayName',TestName,'LineWidth',2)
        hold on

        subplot(2,1,2)
        plot(Time(TimeInd),EffortMea(TimeInd)/MVC*100,'Color',cm(iTest,:),'DisplayName',TestName,'LineWidth',2)
        hold on

        f=figure(2);
        f.Position = [100 300 2000 800];          
        subplot(3,1,1)
        plot(MAVInd,MAV,'o','Color',cm(iTest,:),'DisplayName',TestName)
        hold on
        plot(DroppedInd,DroppedMAV,'*','Color',cm(iTest,:),'LineWidth',2,'DisplayName',TestName)
                ylim([0 0.0005])

        subplot(3,1,2)
        plot(MAVInd,FiltMAV,'Color',cm(iTest,:),'DisplayName',TestName)
        hold on
        plot(DroppedInd,DroppedMAV,'*','Color',cm(iTest,:),'DisplayName',TestName)
        ylim([0 0.0001])

        subplot(3,1,3)
        plot(MAVInd,AmpModulMAV,'Color',cm(iTest,:),'DisplayName',TestName)
        hold on
%         plot(DroppedInd,DroppedMAV/MVC*100,'*','Color',cm(iTest,:),'DisplayName',TestName)

    end
end

figure(1)
subplot(2,1,1)
grid on
legend()
title(ttl)
ylabel("Measured Effort")

subplot(2,1,2)
grid on
legend()
ylabel("Percent MVC")
title(ttl)

figure(2)
subplot(3,1,1)
ylabel("MAV")
legend()
% subplot(3,1,2)
ylabel("BandPass filtered MAV")
legend()
title(ttl)
% ylim([0 0.001])
subplot(3,1,3)
grid on
ylabel("Norm. MAV")
legend()
% ylim([0 .3])


%% Plotting Dropped Frames with Voli only trials
% TestFolders=["jan7" "jan11" "jan12"];
close all
cm = lines(length(TestFolders));
Ref_MVC=0;
% =[10 20 30];
FiltLabel="Unfilt";
ZeroStimFiltLabel="Unfilt";

for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    DroppedsMVC=S.(TestLabel).(ExpLabel).StimMVCVec(2:end);
    vMVC=S.(TestLabel).(ExpLabel).VoliMVCVec;
    sMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('sMVC'));  
    Filts=S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table.('Filt');  

    IndZeroStim=sMVCLevs==Ref_MVC & ZeroStimFiltLabel==Filts ;
    vMVCLevs=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndZeroStim,:).('vMVC'));  
    MAVMean_Reps=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndZeroStim,:).('MAV_Mean'));
    Amp_Modul_Mean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndZeroStim,:).('Amp_Mean'));
    
    sMVCzero_Reps=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndZeroStim,:).('sMVC'));
    vMVCVal_Reps=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndZeroStim,:).('vMVC'));
    
        % Voli Only trials
%     Ind_Reps=(sMVCzero_Reps==Ref_MVC);
    p1=polyfit(vMVCLevs,MAVMean_Reps,1);  %% Combine all the experiments here

        % Dropped frames 
        
    for iStim=1:length(DroppedsMVC)
        IndDropped=sMVCLevs==DroppedsMVC(iStim) & FiltLabel==Filts ;
        DroppedMean=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndDropped,:).('Dropped_Mean'));
        vMVCLevsDropped=(S.(AnaLabel).(ExpLabel).MAV_Mean_reps_table(IndDropped,:).('vMVC'));

        p3=polyfit(vMVCLevsDropped,DroppedMean,1);

        f=figure(1);
        f.Position = [100 300 2000 800];   
        subplot(length(DroppedsMVC),1,iStim)
        plot(vMVCLevs,MAVMean_Reps,'*','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
        hold on
        plot(vMVCLevs,polyval(p1,vMVCLevs),'--','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
        plot(vMVCLevsDropped,DroppedMean,'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',sprintf("%s Dropped",TestFolders(iTest)))
        plot(vMVCLevsDropped,polyval(p3,vMVCLevsDropped),'-','LineWidth',1,'Color',cm(iTest,:),'DisplayName',sprintf("%s Dropped",TestFolders(iTest)))

        ylabel('%% MVC ')
        title(sprintf(' Stim: %d%%', DroppedsMVC(iStim)))

        xlabel('Voli MVC')

    %     ylim([0 0.0001])
    %     subplot(2,1,2)
    %     plot(vMVCVal(Ind),Amp_Modul_Mean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    %     legend
    %     hold on
    %     plot(vMVCVal(Ind),polyval(p2,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    %     ylabel('% Amp modul')
    %     xlabel('% Voli. MVC')
    %     xlim([0 50])
        
    end
    legend
end

%% Plotting fittings 

clear Occ
close all
Log=false;
IndivType=true;
OccType=["Occ_mvc" "Occ_Dropped_mvc" ];
% OccType=["Occ_mvc" "Occ_Dropped_mvc"];
FitVoliPts=[10:60];
Fld=TestFolders(6:end);
for iType=1:length(OccType)
    Type=OccType(iType);

    for iTest=1:length(TestFolders)
        Test=TestFolders(iTest);
        AnaLabel=sprintf("%s_ana",Test);
        EffortType=S.(TestLabel).ExpPar.EffortType;
    lmTable=table();
        
        RowInd=OccTable.('Feat')=="EffortMea" & OccTable.('Filt')=="Unfilt" &...
            OccTable.('Tau')=="3*tau" & OccTable.('MVC_Stim')~="0";
        
        lmTable.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
        lmTable.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
        lmTable.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
        lmTable.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
        lmTable.Test=string((OccTable(RowInd,:).('Test')));
        
        [gVoli, VoliID]=findgroups(lmTable(lmTable.Test==Test,:).('VoliMVC'));
        [gStim, StimID]=findgroups(lmTable(lmTable.Test==Test,:).('StimMVC'));
        TestName=S.(AnaLabel).AnaPar.TestName;

        FitVoliPts=[VoliID(1):VoliID(end)];
        cm=lines(StimID);
        for iStim=1:length(StimID)
            Stim=StimID(iStim);
            StimLabel=sprintf("Stim:%d",Stim);
            
            CoeffRowInd=  OccCoef3.('Log')==string(Log) & OccCoef3.('Indiv')==string(IndivType)...
            & OccCoef3.('Type') == Type & OccCoef3.('Test') == Test;
            
            Coeff_Indiv(1)=OccCoef3(CoeffRowInd,:).('Coeff1');
            Coeff_Indiv(2)=OccCoef3(CoeffRowInd,:).('Coeff2');
            Coeff_Indiv(3)=OccCoef3(CoeffRowInd,:).('Coeff3');
            Coeff_Indiv=double(Coeff_Indiv);

            Occ=zeros(length(VoliID),length(StimID));
            for iVoli=1:length(VoliID)
                OccInd=lmTable.('Test')==Test & lmTable.('StimMVC')==Stim & lmTable.('VoliMVC')==VoliID(iVoli);
                Occ(iVoli,1:sum(OccInd))=lmTable(OccInd,:).('Occ');
            end
            
            FitOcc=Coeff_Indiv(1) + Coeff_Indiv(2)*FitVoliPts + Coeff_Indiv(3)*Stim;
            figure(4)
            subplot(length(TestFolders),length(OccType),length(OccType)*iTest-length(OccType)+iType)
            plot(FitVoliPts,FitOcc,'DisplayName',sprintf("Stim%d",StimID(iStim)),'Color',cm(iStim,:))
            hold on
            plot(VoliID,Occ,'o','DisplayName',sprintf("Stim%d",StimID(iStim)),'Color',cm(iStim,:))
            ylim([-10 50])
            % xlabel(['Percent Voli. Levels'])
            ylabel(TestName)
            % clear Occ
        end
    end
end

% subplot(length(TestFolders),length(OccType),1)
% legend

%% Coefficient plotting

lmTable=table();

iType=2;
Type=OccType(iType);
Log=false;
x_v=10+0.1:0.1:40;
x_s=10+0.1:0.1:40;

for iTest=1:length(TestFolders)
    TestName=TestFolders(iTest);

    EffortType=S.(TestLabel).ExpPar.EffortType;

    RowInd=OccTable.('Feat')=="EffortMea" &...
        OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau" &...
        OccTable.('MVC_Stim')~="0";

    lmTable.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
    lmTable.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
    lmTable.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
    lmTable.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
    lmTable.Test=string((OccTable(RowInd,:).('Test')));

    IndivType=true;
    CoeffRowInd=  OccCoef3.('Log')==string(Log) & OccCoef3.('Indiv')==string(IndivType)...
        & OccCoef3.('Type') == Type & OccCoef3.('Test') == TestName;
    
    Coeff_Indiv(1)=OccCoef3(CoeffRowInd,:).('Coeff1');
    Coeff_Indiv(2)=OccCoef3(CoeffRowInd,:).('Coeff2');
    Coeff_Indiv(3)=OccCoef3(CoeffRowInd,:).('Coeff3');
    
    PVals_Indiv(1)=OccCoef3(CoeffRowInd,:).('pVal1');
    PVals_Indiv(2)=OccCoef3(CoeffRowInd,:).('pVal2');
    PVals_Indiv(3)=OccCoef3(CoeffRowInd,:).('pVal3');
    Coeff_Indiv=double(Coeff_Indiv);
    PVals_Indiv=double(PVals_Indiv);

    IndivType=false;

    CoeffRowInd=  OccCoef3.('Log')==string(Log) & OccCoef3.('Indiv')==string(IndivType)...
        & OccCoef3.('Type') == Type & OccCoef3.('Test') == TestName;
    
    Coeff_Gen(1)=OccCoef3(CoeffRowInd,:).('Coeff1');
    Coeff_Gen(2)=OccCoef3(CoeffRowInd,:).('Coeff2');
    Coeff_Gen(3)=OccCoef3(CoeffRowInd,:).('Coeff3');
    
    PVals_Gen(1)=OccCoef3(CoeffRowInd,:).('pVal1');
    PVals_Gen(2)=OccCoef3(CoeffRowInd,:).('pVal2');
    PVals_Gen(3)=OccCoef3(CoeffRowInd,:).('pVal3');
    Coeff_Gen=double(Coeff_Gen);
    PVals_Gen=double(PVals_Gen);
    
    figure(5)
    subplot(2,2,1)
    plot(1:length(Coeff_Indiv), Coeff_Indiv,'o','DisplayName',TestName,'LineWidth',2)
    hold on
    title('Individual Model Coefficients')

    subplot(2,2,2)
    plot(1:length(PVals_Indiv), PVals_Indiv,'o','DisplayName',TestName,'LineWidth',2)
    hold on
    legend()
    title('Individual Model p-Val')
    
    subplot(2,2,3)
    plot(1:length(Coeff_Gen), Coeff_Gen,'o','DisplayName',TestName,'LineWidth',2)
    hold on
    title('Generalized Model Coefficients')

    subplot(2,2,4)
    plot(1:length(PVals_Gen), PVals_Gen,'o','DisplayName',TestName,'LineWidth',2)
    hold on
    legend()
    title('Generalized Model p-Val')


end

subplot(2,2,2)
plot([1 3],[.05 .05 ],'LineWidth',2,'Color','r','DisplayName','Alpha')
subplot(2,2,4)
plot([1 3],[.05 .05 ],'LineWidth',2,'Color','r','DisplayName','Alpha')


%% Plotting Effort
close all
clc
lbl='Occ';
PlotVoli=2;
PlotStim=4;
% TestFolders=["jan7" "jan11" "jan12"];
TimeRange=[1 12];

cm=lines(6);
ln=["--" ":" "-."];
FiltLabel="GS";
for iType=1:length(OccType)

    TypeLabel=OccType(iType);

    for iTest=1:length(TestFolders)
        AnaLabel=sprintf('%s_ana',TestFolders(iTest));
        TestLabel=sprintf('%s_test',TestFolders(iTest));
        EffortType=S.(TestLabel).ExpPar.EffortType;

        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        FrameInd=stim_freq*TimeRange(1): stim_freq*TimeRange(2);

        ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(lbl);
        RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
        StimMVCVec=S.(TestLabel).(ExpLabel).StimMVCVec;
        VoliMVCVec=S.(TestLabel).(ExpLabel).VoliMVCVec;

        sMVC=StimMVCVec(PlotStim);
        vMVC=VoliMVCVec(PlotVoli);
        IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
        ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);

        RowInd=OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')=="EffortMea" &...
            OccTable.('MVC_Voli')==num2str(vMVC) & OccTable.('MVC_Stim')==num2str(sMVC) &...
            OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau";

        Fv_prime=str2double(OccTable(RowInd,:).('vprime_mvc'));

        E_e_MAV=zeros(length(FrameInd),length(IndTrials));
        E_MAV=zeros(length(FrameInd),length(IndTrials));
        E_MAV_Indv=zeros(length(FrameInd),length(IndTrials));
        E_o_Indv=zeros(length(FrameInd),length(IndTrials));
        E_e_Amp=zeros(length(FrameInd),length(IndTrials));
        E_Amp=zeros(length(FrameInd),length(IndTrials));
        E_Amp_Indv=zeros(length(FrameInd),length(IndTrials));
        E_o=zeros(length(FrameInd),length(IndTrials));
        E_f=zeros(length(FrameInd),length(IndTrials));
        E_fv_prime=zeros(length(FrameInd),length(IndTrials));

        for iTrial=1:length(IndTrials)

            TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
            E_e_MAV(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_e_MAV');
            E_MAV(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_MAV');
            E_MAV_Indv(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_MAV_Indv');

            E_o_Indv(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_o_Indv');
            E_e_Amp(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_e_Amp');
            E_Amp(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Amp');
            E_Amp_Indv(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Amp_Indv');

            E_o(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_o');   
            E_f(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_f');
            E_fv_prime(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Fv_prime');  
        end

        % With averaged (across participants) maintained level
        figure(1)
        subplot(length(OccType),length(TestFolders),iTest+(iType*length(TestFolders)-length(TestFolders)))
        plot(FrameInd,mean(E_f,2),'DisplayName',sprintf('E_f, Trial: %d, %d, %d',IndTrials),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        hold on
        plot(FrameInd,mean(E_e_MAV,2),'DisplayName',sprintf('E_e, Trial: %d, %d, %d',IndTrials),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_MAV,2),'DisplayName',sprintf('E_e+E_o ,Trial: %d, %d, %d',IndTrials),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_fv_prime,2),'DisplayName',sprintf('E_{vprime} ,Trial: %d, %d, %d',IndTrials),'Color',cm(4,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_MAV_Indv,2),'DisplayName',sprintf('E_e+E_{o, indv} ,Trial: %d, %d, %d',IndTrials),'Color',cm(5,:))%,'LineStyle',ln(iTrial))
        title(TestFolders(iTest))
        ylabel(TypeLabel)
        xlabel('Frames')
        grid on

        figure(2)
        subplot(length(OccType),length(TestFolders),iTest+(iType*length(TestFolders)-length(TestFolders)))
        plot(FrameInd,mean(E_f,2),'DisplayName',sprintf('E_f, Trial: %d, %d, %d',IndTrials),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        hold on
        plot(FrameInd,mean(E_e_Amp,2),'DisplayName',sprintf('E_e, Trial: %d, %d, %d',IndTrials),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_Amp,2),'DisplayName',sprintf('E_e+E_o ,Trial: %d, %d, %d',IndTrials),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_fv_prime,2),'DisplayName',sprintf('F_{vprime} ,Trial: %d, %d, %d',IndTrials),'Color',cm(4,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,mean(E_Amp_Indv,2),'DisplayName',sprintf('E_e+E_{o, indv} ,Trial: %d, %d, %d',IndTrials),'Color',cm(5,:))%,'LineStyle',ln(iTrial))
        title(TestFolders(iTest))
        ylabel(TypeLabel)
        xlabel('Frames')
        grid on
    end
    figure(2)
    legend('Location','NorthWest')
end


%%%%%----------------------------------------------

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