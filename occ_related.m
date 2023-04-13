%% Occlusion Analysis
% This file contains all the occlusion related analyses
%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);
%% Plotting the Occ Trials

close all
clc
lbl='Occ';
PlotVoli=2;
PlotStim=4;
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 ];
sMVC=StimMVCLevels(PlotStim);
vMVC=VoliMVCLevels(PlotVoli);

for iTest=1:length(TestFolders)
    AnaStruct=sprintf('%s_ana',TestFolders(iTest));
    TestStruct=sprintf('%s_test',TestFolders(iTest));
    
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(lbl);
    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;

    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

    StimOff=S.(TestStruct).(ExpLabel).StimProfile;
    PlotRange=[ StimOff-10 StimOff+1];
    PreStimOffRange= [StimOff-1 StimOff ];
    PostStimOffRange=[StimOff+.2 StimOff+.3];
    ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);
    
    for iTrial=1:length(IndTrials)

        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        FiltForce=S.(AnaStruct).(ExpLabel).(TrialLabel).data.('Filt_Force');
        Time=S.(AnaStruct).(ExpLabel).(TrialLabel).data.('Time');        
%         
%         PreStimOffInd=PreStimOffRange(1)<=Time & Time<=PreStimOffRange(2);
%         PostStimOffInd=PostStimOffRange(1)<=Time & Time<=PostStimOffRange(2);
%         
        figure(1)
        subplot(length(TestFolders),1,iTest)
        plot(Time,FiltForce,'DisplayName',sprintf('Trial: %d',IndTrials(iTrial)))
        hold on
        
%         PreStimOffForce(iTrial)=mean(FiltForce(PreStimOffInd ));
%         PostStimOffForce(iTrial)=mean(FiltForce(PostStimOffInd ));
    end
    
    title(ttl)
    legend('Location','NorthWest')
    grid on

end

%% Time Constants of RCCurve Trials

%To Do List
% 1- remove filtering of force
% 2- indicing simplifications 
% 3- Saving fixes: dont save the averages, save all the coefficients ^
% 4- Plotting presentations
% 5- make this part of main analysis

% # Filter Design
d1 = designfilt("lowpassiir",'FilterOrder',3, ...
    'HalfPowerFrequency',0.01,'DesignMethod',"butter"); %,'SampleRate',fs

% fvtool(d1)

LPPass = 30;
LPStop = 100;
Ap = .1;
ExpLabel=["RCCurveTrials"];

for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));

    
    fs=S.(TestLabel).ExpPar.fs;  

    d1 = designfilt('lowpassfir','PassbandFrequency',LPPass,...
      'StopbandFrequency',LPStop,'PassbandRipple',Ap,...
      'DesignMethod', 'kaiserwin','SampleRate',fs);
    
    Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign=d1;

end

% # Time Constant Estimation 
% This code will execute the estimation of the time constant 
% Redo Trials must be incorporated 

clc
TauTime=.37;
clear F_mat IndTau
ExpLabel=["RCCurveTrials"];

for iTest=1:length(TestFolders) 
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    
    NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
    DataInd= S.(TestLabel).ExpPar.DataInd;
    DataVars=string(DataInd.Properties.VariableNames);
    iForce=(S.(TestLabel).ExpPar.DataInd.("Force"));
    iTime=(S.(TestLabel).ExpPar.DataInd.("Time"));
    
    PWVal=S.(TestLabel).(ExpLabel).PWTrials;
    TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
    AvgRangePostOff=[TurnOffTime+.5 TurnOffTime+.8];
    AvgRangePreOff=[TurnOffTime-1.5 TurnOffTime];    
    fs=S.(TestLabel).ExpPar.fs;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf("Trial_%d",iTrial);

        F=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Force');
        T=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Time');
         %% 
         %There is no need for this 
         %%
        FiltDesign=Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign;
%         F_filtered = filtfilt(FiltDesign,F);
        F_filtered=F;
        
        Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
        PreAvgForce= mean(F_filtered(Ind));
        Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
        PostAvgForce=mean(F_filtered(Ind));
%         t_end=10.1-TurnOffTime;

        tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1))-PostAvgForce,-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce-PostAvgForce);
%         tau=-t_end/log((PreAvgForce-PostAvgForce)/PreAvgForce);

%         F_tau=TauTime*abs(PreAvgForce-PostAvgForce)+PostAvgForce;
%         Ind=T> TurnOffTime & T<AvgRangePostOff(1); % Ind for Force Drop region 
%         [~,IndTau]=min(abs(F_tau-F_filtered(Ind)));
        
        IndTau=round((TurnOffTime+tau)*fs);
%         IndTau=max((find(abs(PreAvgForce-F_filtered)<=0.63*(PreAvgForce-PostAvgForce))));
%         tau=T(IndTau)-TurnOffTime;
        
        Taus(iTrial,7)= PostAvgForce;
        Taus(iTrial,6)= PreAvgForce;
        Taus(iTrial,5)= iTrial;
        Taus(iTrial,4)= iTest;
%         Taus(iTrial,4)= [TurnOffTime AvgRangePostOff];

        Taus(iTrial,3)= tau;
        Taus(iTrial,2)= IndTau;
        Taus(iTrial,1)= PWVal(iTrial);
        
        F_mat(:,1)=T';
        VarNames(1)={'Time'};
        F_mat(:,iTrial+1)=F_filtered()';
        VarNames{iTrial+1}=char(TrialLabel);

    end
    
    TauVars={'PW','TauInd','TimeConst' ,'Test','TrialNum','PreAvgForce','PostAvgForce'};
    Tau_table=array2table(Taus,'VariableNames',TauVars);
    Redo= false(height(Tau_table),1);
    T_redo = table(Redo,'VariableNames',"Redo");
    Tau_table= [Tau_table T_redo];
    Tau.(TestLabel).(ExpLabel).TimeCons.TurnOffTime=TurnOffTime;
    Tau.(TestLabel).(ExpLabel).TimeCons.AvgRangePostOff=AvgRangePostOff;
    Tau.(TestLabel).(ExpLabel).TimeCons.AvgRangePreOff=AvgRangePreOff;
    Tau.(TestLabel).(ExpLabel).TimeCons.Taus=Taus;
    Tau.(TestLabel).(ExpLabel).TimeCons.Tau_table=Tau_table;
    
    F_table=array2table(F_mat,'VariableNames',VarNames);
    clear F_mat
    Tau.(TestLabel).(ExpLabel).TimeCons.F_filtered=F_table;

    %%
    % Similar algo for the redo trials
    %%
    
    clear F_mat VarNames Taus IndTau

    RedoTrials=S.(TestLabel).(ExpLabel).RedoTrials;
    S.(TestLabel).(ExpLabel).TimeConsRedo.RedoTrials=RedoTrials;
    PWVal=S.(TestLabel).(ExpLabel).PWTrials

    if ~isempty(RedoTrials)
        for iTrial=1:length(RedoTrials)
            TrialLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));

            F=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Force');
            T=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Time');
            
            FiltDesign=Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign;
%             F_filtered = filtfilt(FiltDesign,F);
            F_filtered=F;

            Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
            PreAvgForce= mean(F_filtered(Ind));
            Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
            PostAvgForce=mean(F_filtered(Ind));
            
%             t_end=mean(AvgRangePreOff)-TurnOffTime;
%             tau=-t_end*log((PreAvgForce-PostAvgForce)/PreAvgForce);
            tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1))-PostAvgForce,-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce-PostAvgForce);

%             F_tau=TauTime*abs(PreAvgForce-PostAvgForce)+PostAvgForce;
%             Ind=T> TurnOffTime & T<AvgRangePostOff(1); % Ind for Force Drop region 
%             [~,IndTau]=min(abs(F_tau-F_filtered(Ind)));

            IndTau=round((TurnOffTime+tau)*fs);
        
%             IndTau=IndTau+round(TurnOffTime*fs);
%             tau=T(IndTau)-TurnOffTime;
    %         IndTau=max((find(abs(PreAvgForce-F_filtered)<=0.63*(PreAvgForce-PostAvgForce))));
    %         tau=T(IndTau)-TurnOffTime;
        
            Taus(iTrial,7)= PostAvgForce;
            Taus(iTrial,6)= PreAvgForce;
            Taus(iTrial,5)= RedoTrials(iTrial);
            Taus(iTrial,4)= iTest;
            Taus(iTrial,3)= tau;
            Taus(iTrial,2)= IndTau;
            Taus(iTrial,1)= PWVal(RedoTrials(iTrial));

            F_mat(:,iTrial)=F_filtered()';
            VarNames{iTrial}=char(TrialLabel);
        end

        Tau_table_redo=array2table(Taus,'VariableNames',TauVars);
        
        Redo= true(height(Tau_table_redo),1);
        T_redo = table(Redo,'VariableNames',"Redo");
        Tau_table_redo= [Tau_table_redo T_redo];
        Tau.(TestLabel).(ExpLabel).TimeConsRedo.TurnOffTime=TurnOffTime;
        Tau.(TestLabel).(ExpLabel).TimeConsRedo.AvgRangePostOff=AvgRangePostOff;
        Tau.(TestLabel).(ExpLabel).TimeConsRedo.AvgRangePreOff=AvgRangePreOff;
        Tau.(TestLabel).(ExpLabel).TimeConsRedo.Taus=Taus;
        Tau.(TestLabel).(ExpLabel).TimeConsRedo.Tau_table=Tau_table_redo;

        F_table_redo=array2table(F_mat,'VariableNames',VarNames);
        clear F_mat

        Tau.(TestLabel).(ExpLabel).TimeConsRedo.F_filtered=F_table_redo;
        F_table=[F_table F_table_redo];
        Tau_table=[Tau_table; Tau_table_redo];        

    end
    Tau.(TestLabel).(ExpLabel).Tau_table=Tau_table;
    Tau.(TestLabel).(ExpLabel).F_table=F_table;

end


% # Describe with Stats
% Mean, std
 
clear MeanTau StdTau
% Incorporate redos

for iTest= 1:length(TestFolders)
    TestLabel= sprintf("%s_test",TestFolders{iTest});
    
    Tau_table=Tau.(TestLabel).(ExpLabel).Tau_table;
    RedoInd=table2array(Tau_table(:,"Redo"));
    Redo_table=Tau_table(RedoInd,"TrialNum");
    Tau_table(table2array(Redo_table),:)=Tau_table(table2array(Tau_table(:,"Redo")),:);
    Tau_table(RedoInd,:)=[];
    Tau.(TestLabel).(ExpLabel).Tau_incorp=Tau_table;
end

Tau_stats=table([],[],[],[],'VariableNames',["PW","mean_TimeConst","std_TimeConst","Test"]);
% Finding groups
for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    Tau_{iTest}=Tau.(TestLabel).(ExpLabel).Tau_incorp;
    [ID_PW{iTest},PWPoints{iTest}]=findgroups(Tau_{iTest}(:,1));
%     PWPoints{iTest} = renamevars(PWPoints{iTest},["PW"],[sprintf('%s PW',TestFolders{iTest})]);
    PWPoints{iTest}
end

PWInd{iTest}=ID_PW{iTest} == Ind(iTest);

for iTest=1:length(TestFolders)
    TestLabel=sprintf("%s",TestFolders{iTest});

    for PWPointsInd=1:height(PWPoints{iTest})

        if isempty( PWPointsInd)
           disp('PW value was not found') 
           return
        end

        PWInd{iTest}=ID_PW{iTest} == PWPointsInd;
        MeanTau{iTest}(PWPointsInd,:) = varfun(@mean, Tau_{iTest}(PWInd{iTest},'TimeConst'), 'InputVariables', @isnumeric);
        StdTau{iTest}(PWPointsInd,:) = varfun(@std, Tau_{iTest}(PWInd{iTest},'TimeConst'), 'InputVariables', @isnumeric);

    end
    Tau_stats=[Tau_stats; PWPoints{iTest} MeanTau{iTest} StdTau{iTest}...
        table([TestLabel+strings(height(PWPoints{iTest}),1)],'VariableNames',"Test")];
end

MTau = array2table(ones(1,length(TauVars)+1),'VariableNames',[TauVars "Redo"]);

for iTest=1:length(TestFolders)
    
    MTau=[MTau; Tau_{iTest}];
end
 
MTau(1,:)=[];


% # Change Variable Types of the Tau Table
    MTau = convertvars(MTau,{'Test'},'categorical');
    MTau.('Test')= renamecats(MTau.('Test'),TestFolders)
    Tau_stats
    
% # Export Add to Original File
% Export as a mat file, and table

FileNameExtension="_tau";
for iTest=1:length(TestFolders)
    TestFile=sprintf('%s_test',TestFolders(iTest));
    ChrTest=char(TestFile);
    TestLabel=sprintf("%s",ChrTest(1:end-5));

    FoldLabel=TestFolders{iTest};
    DirLabel=sprintf('%s/%s%s_test',FoldLabel,FoldLabel,FileNameExtension);
    save(DirLabel,'-struct','Tau',string(ChrTest))
    %%
    %Tau table ".csv" save
    DirLabelCSV=sprintf('%s/%s_test%s.csv',FoldLabel,FoldLabel,FileNameExtension);
    writetable( Tau.(ChrTest).(ExpLabel).Tau_incorp, DirLabelCSV)
    %%
    %Tau stats ".csv" save
    DirLabelCSV=sprintf('%s/%s_taustats.csv',FoldLabel,FoldLabel);
    writetable( Tau_stats(Tau_stats.Test==FoldLabel,:), DirLabelCSV)

    
end
%all the results in one file
writetable( MTau,'tau_estimates.csv')

%%
% # ---------------LoaD-------------------
clc
clear all
% TestFiles=["dec5_tau_test","nov28_2_tau_test","nov27_tau_test","nov8_tau_test"];
% TestFolder=["dec5","nov28_2","nov27","nov8"];
TestFiles=["jan7_tau_test" "jan11_tau_test" "jan12_tau_test"];
TestFolders=["jan7" "jan11" "jan12"];
ExpLabel=["RCCurveTrials"];
K = load_test(TestFolders,TestFiles);

for iTest=1:length(TestFiles)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    Tau{iTest}=K.(TestLabel).(ExpLabel).Tau_table;
    [ID_PW{iTest},PWPoints{iTest}]=findgroups(Tau{iTest}(:,1));
    PWPoints{iTest} = renamevars(PWPoints{iTest},["PW"],[sprintf('%s PW',TestFolders{iTest})]);
    PWPoints{iTest}
    
    %load tables 
    DirLabelCSV=sprintf('%s/%s_taustats.csv',TestFolders{iTest},TestFolders{iTest});
    Tau_Cell{iTest}= readtable(DirLabelCSV);

end

DirLabelCSV=['tau_estimates.csv'];

Tau_table= readtable(DirLabelCSV);

	% load the main analysis  
for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders(iTest));
end

S = load_test(TestFolders,TestFiles);

%% # Plotting 
clc
PWPoints{1:length(TestFolders)}
close all
PlotPW = [130 120 120 100 100 105];
for iTest=1:length(TestFolders)
    TestLabel=sprintf('%s_test',TestFolders{iTest});

    Ind{iTest}=find(table2array(PWPoints{iTest})==PlotPW(iTest));

    if isempty( Ind{iTest})
        sprintf('PW value was not found for %s',TestFolders{iTest})
        return
    end

    PWInd{iTest}=ID_PW{iTest} == Ind{iTest};
    Tau{iTest}(PWInd{iTest},:);
    TrialNums{iTest}=Tau{iTest}(:,'TrialNum');

    PWInd_F=logical([1 PWInd{iTest}']');
    F_PW=K.(TestLabel).(ExpLabel).F_table(:,PWInd_F);
    taus=Tau_table.('TimeConst')(PWInd{iTest});
    trialnum=Tau_table.('TrialNum')(PWInd{iTest});

    %%
    %Plotting here
    h= figure(1);
    set(h, 'Visible', 'on');
    subplot(length(TestFolders),1,iTest)
    plot(table2array(F_PW(:,1)),table2array(F_PW(:,2:end)),'LineWidth',2)
    hold on
    plot([0 5 6 10 10.00001 11],[0 0 PlotPW(iTest) PlotPW(iTest) 0  0]/20)
    xlabel('Time(s)')
    ylabel('Force(N)')
    PWvalues=table2array(PWPoints{iTest});
    TrialVals=table2array(TrialNums{iTest}(PWInd{iTest},'TrialNum'));
    trials_str1=sprintf("Trials: " );
    trials_str2=sprintf("%d, ",TrialVals );
    PW_str=sprintf("PW: %d ",PlotPW(iTest) );
    str_test=sprintf("Test: %s ",TestFolders{iTest} );
    title(strcat(trials_str1,trials_str2,PW_str,",",str_test))
    lgd_str{1}= sprintf ("Trial %d tau=%.2f",TrialVals(1),taus(1));
    lgd_str{2}= sprintf ("Trial %d tau=%.2f",TrialVals(2),taus(2));
    lgd_str{3}= sprintf ("Trial %d tau=%.2f",TrialVals(3),taus(3));
    lgd_str{4}= sprintf ("Norm PW");

    legend(lgd_str,'Location','NorthEast','AutoUpdate','off')
    xlim([5 11.9])
    grid on
    

    %% 
    %annotations and mark important points
    %%
    a = get(gca,'Children');
    y1data = get(a, 'YData');
    y1min=min( [min(y1data{1}) min(y1data{2}) min(y1data{3}) ]);
    y1max=max( [max(y1data{1}) max(y1data{2}) max(y1data{3}) ]);

    PreAvgForce=table2array(Tau{iTest}(PWInd{iTest},'PreAvgForce'));
    PostAvgForce=table2array(Tau{iTest}(PWInd{iTest},'PostAvgForce'));
    AvgRangePreOff=K.(TestLabel).(ExpLabel).TimeCons.AvgRangePreOff;
    AvgRangePostOff=K.(TestLabel).(ExpLabel).TimeCons.AvgRangePostOff;
    TurnOffTime=K.(TestLabel).(ExpLabel).TimeCons.TurnOffTime;

    plot([TurnOffTime TurnOffTime]',[0 y1max]','--','Color','k')
    plot([AvgRangePreOff(1) AvgRangePreOff(2)],[PreAvgForce PreAvgForce],'--','Color','k')
    plot([AvgRangePostOff(1) AvgRangePostOff(2)],[PostAvgForce PostAvgForce],'--','Color','k')

    

end
%%
% # Determining Force, MAV Occlusion 
% Occ eqn F_o= Fs-Fs'   
% Fs : RC curve value with the corresponding PW -> R(PW)
% Fs': Force value at 3*tau

% To do :
% 1- incorporate redos: already incorporated earlier 
% 2- a section goes to main analysis
% 3- what to do with mav target

tau = mean(Tau_table.('TimeConst')); % average time constant to be used 
OccRefs = [3*tau 4*tau 5*tau]; % referance time for occlusion to be calculated
TauLabels=["3*tau" "4*tau" "5*tau"];

AnaLabel=sprintf("%s_ana",TestFolders(1));
 % Fitting the RC curve: turns out no need for this 
lbl=S.(AnaLabel).AnaPar.ExpTable.('RC');
PWPoints=S.(TestLabel).(lbl).PWPoints;
RCVar=S.(TestLabel).(lbl).RCVar;
MVClevels= (PWPoints-min(PWPoints))/max(PWPoints-min(PWPoints))*30;
MeanForce=S.(TestLabel).(lbl).MeanForce;

Occ=[ [] [] [] [] [] [] [] [] [] [] [] []];

for iTau=1:length(OccRefs)
    OccRef=OccRefs(iTau);

    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders(iTest));
        AnaLabel=sprintf("%s_ana",TestFolders(iTest));
        FiltLabels=S.(AnaLabel).AnaPar.FiltLabels;
        
        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        OccRefMargin= 1/stim_freq ; % in secs
        PreOffMargin=10/stim_freq; % Avg period before stim turning off 
    %----This section should go to main_analysis 
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('MVC');
        FiltLabel="Unfilt";
        MVC=S.(TestLabel).(ExpLabel).MVC;
        for iTrial=1:length(S.(TestLabel).(ExpLabel).NumofTrials)
            TrialLabel=sprintf('Trial_%d',iTrial );
            TrialsMAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(25:65,:).('MAV_vEMG'));
        end
        
        MAV_max=mean(TrialsMAV);
    %----
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('Occ');
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat,...
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);
        
%         ForceOcc=[ [] [] [] [] [] [] [] []];
        for iTrial=1:NumofTrials
            TrialLabel=sprintf("Trial_%d", iTrial);
            
            Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');
            PostOffInd= TurnOffTime+OccRef-OccRefMargin<=Time & TurnOffTime+OccRef+OccRefMargin>=Time;
            PreOffInd=TurnOffTime-PreOffMargin<=Time & TurnOffTime>=Time;
            
            Force=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Filt_Force');

            [~,TurnOffInd]=min(abs(Time-TurnOffTime));
            PostForceLevel= mean(Force( PostOffInd));
            PreForceLevel= mean(Force( PreOffInd));
            PW_level=RepMatTable(iTrial,:).('PW');
            F_v=RepMatTable(iTrial,:).('Voli_Force');
            F_stim=RepMatTable(iTrial,:).('Stim_Force');
            F_target=RepMatTable(iTrial,:).('Target_Level');
    
%             F_s=gompertz(RCVar,PW_level);

            F_vprime=PostForceLevel;
            F_occ=F_vprime-F_v;
            
            F_occ_mvc=F_occ/MVC*100;
            F_vprime_mvc=F_vprime/MVC*100;
            F_v_mvc=F_v/MVC*100;
            F_target_mvc=F_target/MVC*100;

            Occ=[ Occ;  F_occ  F_vprime F_v F_target F_occ_mvc F_vprime_mvc...
                F_v_mvc F_target_mvc TestFolders(iTest) TauLabels(iTau) "Force" "Unfilt" table2array(RepMatTable(iTrial,:) )];
        
            PostOffFrameInd=round((TurnOffTime+OccRef-OccRefMargin)*stim_freq:( TurnOffTime+OccRef+OccRefMargin)*stim_freq);
            PreOffFrameInd=round((TurnOffTime-PreOffMargin)*stim_freq:( TurnOffTime)*stim_freq);
            MVC_Voli=RepMatTable(iTrial,:).('MVC_Voli');
            MVC_Stim=RepMatTable(iTrial,:).('MVC_Stim');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('vMVC');
            MAV_v=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps(MVC_Voli==vMVC & sMVC==0,:).('Mean');
            
            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels{iFilt};
                
                MAV_target=F_target; %%% ---- Fix This 

                PreOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(PreOffFrameInd,:).('MAV_vEMG'));
                PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(PostOffFrameInd,:).('MAV_vEMG'));
                MAV_vprime=PostOffMAV;
                MAV_occ=MAV_vprime-MAV_v;
                MAV_occ_mvc=MAV_occ/MAV_max*100;
                MAV_vprime_mvc=MAV_vprime/MAV_max*100;
                MAV_v_mvc=MAV_v/MAV_max*100;
                MAV_target_mvc=MAV_target/MAV_max*100;
                
                Occ=[ Occ;  MAV_occ  MAV_vprime MAV_v MAV_target MAV_occ_mvc ...
                    MAV_vprime_mvc MAV_v_mvc MAV_target_mvc TestFolders(iTest) TauLabels(iTau) "MAV" FiltLabel table2array(RepMatTable(iTrial,:)) ];

            end
        end
    end
end

OccTable=table(Occ(:,1),Occ(:,2),Occ(:,3),Occ(:,4),Occ(:,5), Occ(:,6), Occ(:,7),...
    Occ(:,8),Occ(:,9),Occ(:,10),Occ(:,11),Occ(:,12),Occ(:,13),Occ(:,14),Occ(:,15),Occ(:,16),Occ(:,17),Occ(:,18),'VariableNames',...
    [ "Occ" "vprime" "v" "Target" "Occ_mvc" "vprime_mvc" "v_mvc" "Target_mvc" "Test" "Tau" "Feat" "Filt" "Target_Level" " Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW"])

writetable( OccTable,'occlusion_v3.csv')


%% Occlusion from dropped frames 

DirLabelCSV=sprintf('%s/%s_taustats.csv',TestFolders{iTest},TestFolders{iTest});
Tau_Cell{iTest}= readtable(DirLabelCSV);

DirLabelCSV=['tau_estimates.csv'];





