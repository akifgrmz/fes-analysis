%% Occlusion Analysis
% This file contains all the occlusion related analyses
%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "apr20" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% Plotting the Occ Trials

close all
clc
lbl='Occ';
PlotVoli=3;
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
% 1- remove filtering of force^
% 2- indicing simplifications ^
% 3- Saving fixes: dont save the averages, save all the coefficients ^
% 4- Plotting presentations^
% 5- make this part of main analysis
% 6- EMG t
% # Filter Design
d1 = designfilt("lowpassiir",'FilterOrder',3, ...
    'HalfPowerFrequency',0.01,'DesignMethod',"butter"); %,'SampleRate',fs

% fvtool(d1)

LPPass = 30;
LPStop = 100;
Ap = .1;
ExpLabel=["RCCurveTrials"];

% for iTest=1:length(TestFolders)
%     TestLabel=sprintf('%s_test',TestFolders(iTest));
%     AnaLabel=sprintf('%s_ana',TestFolders(iTest));
% 
%     fs=S.(TestLabel).ExpPar.fs;  
% 
%     d1 = designfilt('lowpassfir','PassbandFrequency',LPPass,...
%       'StopbandFrequency',LPStop,'PassbandRipple',Ap,...
%       'DesignMethod', 'kaiserwin','SampleRate',fs);
%     
%     Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign=d1;
% 
% end

% # Time Constant Estimation 
% This code will execute the estimation of the time constant 
% Redo Trials must be incorporated 

clc
TauTime=.37;
clear F_mat IndTau
ExpLabel=["RCCurveTrials"];
Featlabels=["Force" ];
% for iFeat =1:length(FeatLabels)
%     FeatLabel=Featlabels(iFeat);
    
for iTest=1:length(TestFolders) 
    TestLabel=sprintf("%s_test",TestFolders(iTest));
    AnaLabel=sprintf("%s_ana",TestFolders(iTest));

    NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
    DataInd= S.(TestLabel).ExpPar.DataInd;
    DataVars=string(DataInd.Properties.VariableNames);
    iForce=(S.(TestLabel).ExpPar.DataInd.("Force"));
    iTime=(S.(TestLabel).ExpPar.DataInd.("Time"));
    
    PWVal=S.(TestLabel).(ExpLabel).PWTrials;
    TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
    AvgRangePostOff=[TurnOffTime+.5 TurnOffTime+.8];
    AvgRangePreOff=[TurnOffTime-1.5 TurnOffTime]; 
    
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    PostOffFrameRange=round([AvgRangePostOff(1)*stim_freq AvgRangePostOff(2)*stim_freq]);
    PostOffFrameRangeInd=[PostOffFrameRange(1): PostOffFrameRange(2)];
    PreOffFrameRange=round([AvgRangePreOff(1)*stim_freq AvgRangePreOff(2)*stim_freq]);
    PreOffFrameRangeInd=[PreOffFrameRange(1): PreOffFrameRange(2)];
    fs=S.(TestLabel).ExpPar.fs;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf("Trial_%d",iTrial);

        F=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Force');
        T=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Time');
         %% 
         %There is no need for this 
         %%
%         FiltDesign=Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign;
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
    
%     clear F_mat VarNames Taus IndTau
% 
%     RedoTrials=S.(TestLabel).(ExpLabel).RedoTrials;
%     S.(TestLabel).(ExpLabel).TimeConsRedo.RedoTrials=RedoTrials;
%     PWVal=S.(TestLabel).(ExpLabel).PWTrials
% 
%     if ~isempty(RedoTrials)
%         for iTrial=1:length(RedoTrials)
%             TrialLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));
% 
%             F=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Force');
%             T=S.(TestLabel).(ExpLabel).(TrialLabel).data.('Time');
%             
%             FiltDesign=Tau.(TestLabel).(ExpLabel).TimeCons.FiltDesign;
% %             F_filtered = filtfilt(FiltDesign,F);
%             F_filtered=F;
% 
%             Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
%             PreAvgForce= mean(F_filtered(Ind));
%             Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
%             PostAvgForce=mean(F_filtered(Ind));
%             
% %             t_end=mean(AvgRangePreOff)-TurnOffTime;
% %             tau=-t_end*log((PreAvgForce-PostAvgForce)/PreAvgForce);
%             tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1))-PostAvgForce,-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce-PostAvgForce);
% 
% %             F_tau=TauTime*abs(PreAvgForce-PostAvgForce)+PostAvgForce;
% %             Ind=T> TurnOffTime & T<AvgRangePostOff(1); % Ind for Force Drop region 
% %             [~,IndTau]=min(abs(F_tau-F_filtered(Ind)));
% 
%             IndTau=round((TurnOffTime+tau)*fs);
%         
% %             IndTau=IndTau+round(TurnOffTime*fs);
% %             tau=T(IndTau)-TurnOffTime;
%     %         IndTau=max((find(abs(PreAvgForce-F_filtered)<=0.63*(PreAvgForce-PostAvgForce))));
%     %         tau=T(IndTau)-TurnOffTime;
%         
%             Taus(iTrial,7)= PostAvgForce;
%             Taus(iTrial,6)= PreAvgForce;
%             Taus(iTrial,5)= RedoTrials(iTrial);
%             Taus(iTrial,4)= iTest;
%             Taus(iTrial,3)= tau;
%             Taus(iTrial,2)= IndTau;
%             Taus(iTrial,1)= PWVal(RedoTrials(iTrial));
% 
%             F_mat(:,iTrial)=F_filtered()';
%             VarNames{iTrial}=char(TrialLabel);
%         end
% 
%         Tau_table_redo=array2table(Taus,'VariableNames',TauVars);
%         
%         Redo= true(height(Tau_table_redo),1);
%         T_redo = table(Redo,'VariableNames',"Redo");
%         Tau_table_redo= [Tau_table_redo T_redo];
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.TurnOffTime=TurnOffTime;
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.AvgRangePostOff=AvgRangePostOff;
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.AvgRangePreOff=AvgRangePreOff;
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.Taus=Taus;
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.Tau_table=Tau_table_redo;
% 
%         F_table_redo=array2table(F_mat,'VariableNames',VarNames);
%         clear F_mat
% 
%         Tau.(TestLabel).(ExpLabel).TimeConsRedo.F_filtered=F_table_redo;
%         F_table=[F_table F_table_redo];
%         Tau_table=[Tau_table; Tau_table_redo];        
% 
%     end

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
writetable( MTau,'tau_estimates2.csv')

%%
% # ---------------LoaD-------------------
clc
clear all
% TestFiles=["dec5_tau_test","nov28_2_tau_test","nov27_tau_test","nov8_tau_test"];
% TestFolder=["dec5","nov28_2","nov27","nov8"];

TestFolders=["jan7" "jan11" "jan12" "apr20" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_tau_test",TestFolders{iTest});
end
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

DirLabelCSV=['tau_estimates2.csv'];

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
%% Occlusion Analysis
% # Determining Force, MAV Occlusion 
% Occ eqn F_o= Fs-Fs'   
% Fs : RC curve value with the corresponding PW -> R(PW)
% Fs': Force value at 3*tau

% To do :
% 1- incorporate redos: already incorporated earlier 
% 2- a section goes to main analysis
% 3- what to do with mav target
% 4- Time const. for MAV drop

tau = mean(Tau_table.('TimeConst')); % average time constant to be used 
OccRefs = [3*tau 4*tau 5*tau]; % referance time for occlusion to be calculated
TauLabels=["3*tau" "4*tau" "5*tau"];
% TestFolders=["jan7" "jan11" "jan12"];% "apr20"];


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
        FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
        
        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        OccRefMargin= 1/stim_freq ; % in secs
        PreOffMargin=10/stim_freq; % Avg period before stim turning off 
    %----This section should go to main_analysis 
    % MAV_MAX = avg MAV at the last 2 seconds
        AvgTime=2;
        AvgInd=stim_freq*AvgTime-5;
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('MVC');
        FiltLabel="Unfilt";
        MVC=S.(TestLabel).(ExpLabel).MVC;
        for iTrial=1:S.(TestLabel).(ExpLabel).NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial );
            TrialsMAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(end-AvgInd:end,:).('MAV_vEMG'));
            TrialsAmp_MAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(end-AvgInd:end,:).('Amp_MAV_vEMG'));

        end
        MAV_max=mean(TrialsMAV);
        S.(AnaLabel).(ExpLabel).MAV_MAX=MAV_max;
        
        AmpModul_MAV_max=mean(TrialsAmp_MAV);
        S.(AnaLabel).(ExpLabel).AmpModul_MAV_MAX=AmpModul_MAV_max;

        % Theoretical MAV_MAX
        sMVC=0;
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('Occ');
        sMVCzero_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
        sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC'); 
        vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');

        MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean.('Mean');
        Amp_Modul_Mean=S.(AnaLabel).(ExpLabel).Amp_Modul_Mean.('Mean');

        Ind=(sMVCzero==sMVC);
        Ind_Reps=(sMVCzero_Reps==sMVC);
        
        poly1=polyfit(vMVCVal(Ind),MAVMean(Ind),1);
        poly2=polyfit(vMVCVal(Ind),Amp_Modul_Mean(Ind),1);

        MAV_max_theo=polyval(poly1,100);
        AmpModul_max_theo=polyval(poly2,100);

        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('MVC');
        S.(AnaLabel).(ExpLabel).MAV_MAX_theo=MAV_max_theo;
        S.(AnaLabel).(ExpLabel).AmpModul_max_theo=AmpModul_max_theo;

        AvgTime=1.5;
        AvgInd=round(stim_freq*AvgTime:stim_freq*(AvgTime+3));
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('Occ');
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
    %----
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('Occ');
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat,...
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);
        
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
            AmpModul_MAV_v=S.(AnaLabel).(ExpLabel).Amp_Modul_Mean_Reps(MVC_Voli==vMVC & sMVC==0,:).('Mean');

            
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
            
                AmpModul_MAV_target=F_target; %%% ---- Fix This 
                
                AmpModul_PreOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(PreOffFrameInd,:).('Amp_MAV_vEMG'));
                AmpModul_PostOffMAV=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(PostOffFrameInd,:).('Amp_MAV_vEMG'));
                AmpModul_MAV_vprime=AmpModul_PostOffMAV;
                AmpModul_MAV_occ=AmpModul_MAV_vprime-AmpModul_MAV_v;
                AmpModul_MAV_occ_mvc=AmpModul_MAV_occ/AmpModul_MAV_max*100;
                AmpModul_MAV_vprime_mvc=AmpModul_MAV_vprime/AmpModul_MAV_max*100;
                AmpModul_MAV_v_mvc=AmpModul_MAV_v/AmpModul_MAV_max*100;
                AmpModul_MAV_target_mvc=AmpModul_MAV_target/AmpModul_MAV_max*100;
                
                Occ=[ Occ;  AmpModul_MAV_occ  AmpModul_MAV_vprime AmpModul_MAV_v AmpModul_MAV_target AmpModul_MAV_occ_mvc ...
                    AmpModul_MAV_vprime_mvc AmpModul_MAV_v_mvc AmpModul_MAV_target_mvc TestFolders(iTest) TauLabels(iTau) "Amp_Modul" FiltLabel table2array(RepMatTable(iTrial,:)) ];
            
            end
        end
        
        OccTest=Occ(Occ(:,9)==TestFolders(iTest),:);
        OccTest=table(OccTest(:,1),OccTest(:,2),OccTest(:,3),OccTest(:,4),OccTest(:,5), OccTest(:,6), OccTest(:,7),...
            OccTest(:,8),OccTest(:,9),OccTest(:,10),OccTest(:,11),OccTest(:,12),OccTest(:,13),OccTest(:,14),OccTest(:,15),OccTest(:,16),OccTest(:,17),OccTest(:,18),'VariableNames',...
            [ "Occ" "vprime" "v" "Target" "Occ_mvc" "vprime_mvc" "v_mvc" "Target_mvc" "Test" "Tau" "Feat" "Filt" "Target_Level" " Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW"]);
        S.(AnaLabel).(ExpLabel).OccTest=OccTest;
        
        clear OccTest
    end
end

OccTable=table(Occ(:,1),Occ(:,2),Occ(:,3),Occ(:,4),Occ(:,5), Occ(:,6), Occ(:,7),...
    Occ(:,8),Occ(:,9),Occ(:,10),Occ(:,11),Occ(:,12),Occ(:,13),Occ(:,14),Occ(:,15),Occ(:,16),Occ(:,17),Occ(:,18),'VariableNames',...
    [ "Occ" "vprime" "v" "Target" "Occ_mvc" "vprime_mvc" "v_mvc" "Target_mvc" "Test" "Tau" "Feat" "Filt" "Target_Level" " Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW"])

writetable( OccTable,'occlusion_v4.csv')

%% OccMap Interpolation 
clc
clear Occ_MVC OccMap
% occ map: rows: %stim, columns: %voli
% TestFolders=["jan7" "jan11" "jan12" "apr20"];
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    OccTest=S.(AnaLabel).(ExpLabel).OccTest; 
    [~,VoliLevels]=findgroups(OccTest.('MVC_Voli'));
    [~,StimLevels]=findgroups(OccTest.('MVC_Stim'));

    for iVoli=1:length(VoliLevels)
        VoliLabel=VoliLevels(iVoli);
        
        for iStim=2:length(StimLevels)
            StimLabel=StimLevels(iStim);
            
            % 1- Force
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Force'...
                &  OccTest.('Tau')=='3*tau';
            
            OccForce(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));
            OccMap(iTest,iVoli,iStim-1)=OccForce(iVoli,iStim-1);

            % 2- MAV
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='MAV'...
                &  OccTest.('Tau')=='3*tau';
            
            OccMAV(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));
            
            % 3- Amplitude Modul.
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Amp_Modul'...
                &  OccTest.('Tau')=='3*tau';
            
            OccAmpModul(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));

        end
    end
    S.(AnaLabel).(ExpLabel).OccForce=OccForce;
    S.(AnaLabel).(ExpLabel).OccMAV=OccMAV;
    S.(AnaLabel).(ExpLabel).OccAmpModul=OccAmpModul;

    clear OccForce OccMAV OccAmpModul
end

dx=1;
Occ_MVC=squeeze(mean(OccMap));
VoliLevels=str2double(VoliLevels);
StimLevels=str2double(StimLevels(2:end));
vx=[VoliLevels(1):1/(StimLevels(end)-StimLevels(1))/dx:VoliLevels(end)];
vy=[StimLevels(1):1/(VoliLevels(end)-VoliLevels(1))/dx:StimLevels(end)];

% vq=interpn(VoliLevels,StimLevels,Occ_MVC,vx',vy);

% figure
% mesh(vx,vy,vq)
% xlabel('Stim MVC(%)')
% ylabel('Voli MVC(%)')
% zlabel('Occlusion MVC (%)')

%% Effort simulation 
ExpstoAna= ["Occ"];
AnaLabel=sprintf('%s_ana',TestFolders(1));

for iExp=1:length(ExpstoAna)
    ExpLabels(iExp)=S.(AnaLabel).AnaPar.ExpTable.(ExpstoAna(iExp));
end

FeattoAna=["Filt_MAV_vEMG"];
RampTime=[5 10];
ConstTime=[10 15];
TurnOffTime=[15 17];
RampInd=[stim_freq*RampTime(1):stim_freq*RampTime(2)];
ConstInd=[stim_freq*ConstTime(1):stim_freq*ConstTime(2)];
clear Effort

for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf("%s_test",TestFolders(iTest));

    MAV_MAX=S.(AnaLabel).MVCTrials.MAV_MAX_theo;
    AmpModul_max_theo=S.(AnaLabel).MVCTrials.AmpModul_max_theo;
    
    F_MAX=S.(TestLabel).MVCTrials.MVC;

    FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
    OccTest=S.(AnaLabel).(ExpLabel).OccTest; 
    RCLabel=S.(AnaLabel).AnaPar.ExpTable.('RC');
    RCVar=S.(TestLabel).(RCLabel).RCVar; 
    MVCLabel=S.(AnaLabel).AnaPar.ExpTable.('MVC');
    MVC=S.(TestLabel).(MVCLabel).MVC; 

    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels(iExp);

        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;

        for iTrial=1:NumofTrials 
            TrialLabel=sprintf("Trial_%d",iTrial);

            DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
            PWofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames;
            StimMVC=S.(TestLabel).(ExpLabel).RepTableMat(iTrial,5);
            StimMVCofFrames=gompertz(RCVar,PWofFrames)/MVC*100+0.00001;

            for iFilt=1:length(FiltLabels)
                FiltLabel=FiltLabels(iFilt);

                KeepInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd;
                % 1- Force 
                Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames(:,KeepInd))';
                Effort_f=Force/F_MAX*100;

                % 2-MAV
                Filt_Feats=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Filt_MAV_vEMG');
                MAV_Noise=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MAV_Noise;

                Effort_e_MAV=(Filt_Feats-MAV_Noise)/(MAV_MAX-MAV_Noise)*100;
                OccStimInd=round(StimMVCofFrames(KeepInd)*length(vy)/max(StimLevels));
                OccVoliInd=round([Effort_e_MAV]*length(vx)/max(VoliLevels));
                OccVoliInd(OccVoliInd>=length(vx))=length(vx);
                OccVoliInd(OccVoliInd<=0)=1;
                OccStimInd(OccStimInd<=0)=1;
                OccStimInd(OccStimInd>=length(vy))=length(vy);

                for iFrame=1:length(Filt_Feats)
                    Effort_o_MAV(iFrame,1)=vq(OccVoliInd(iFrame),OccStimInd(iFrame));
                end
                
                % 3-Amp Modul
                AmpModul_Feats=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Amp_MAV_vEMG');
                AmpModul_MAV_Noise=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModul_MAV_Noise;

                Effort_e_Amp=(AmpModul_Feats-AmpModul_MAV_Noise)/(AmpModul_max_theo-AmpModul_MAV_Noise)*100;
                OccStimInd=round(StimMVCofFrames(KeepInd)*length(vy)/max(StimLevels));
                OccVoliInd=round([Effort_e_Amp]*length(vx)/max(VoliLevels));
                OccVoliInd(OccVoliInd>=length(vx))=length(vx);
                OccVoliInd(OccVoliInd<=0)=1;
                OccStimInd(OccStimInd<=0)=1;
                OccStimInd(OccStimInd>=length(vy))=length(vy);

                for iFrame=1:length(AmpModul_Feats)
                    Effort_o_Amp(iFrame,1)=vq(OccVoliInd(iFrame),OccStimInd(iFrame));
                end
                
                % 4- Saving the Results

                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Effort_e_MAV')=Effort_e_MAV;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Effort_o_MAV')=Effort_o_MAV;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Effort_f')=Effort_f;

                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Effort_e_Amp')=Effort_e_Amp;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Effort_o_Amp')=Effort_o_Amp;
                S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Effort_f')=Effort_f;

            end
        end
    end
end

%% Plotting Effort
close all
clc
lbl='Occ';
PlotVoli=4;
PlotStim=4;
% TestFolders=["jan7" "jan12" "apr20"];
TimeRange=[1 15];
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 ];
sMVC=StimMVCLevels(PlotStim);
vMVC=VoliMVCLevels(PlotVoli);
cm=lines(3);
ln=["--" ":" "-."];
FiltLabel="Comb";
for iTest=1:4
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=[stim_freq*TimeRange(1): stim_freq*TimeRange(2)];
    
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(lbl);
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
    ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);
    
    for iTrial=1:length(IndTrials)
    
        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        Effort_e=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).('Effort_e_MAV');
        Effort_o=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).('Effort_o_MAV');   
        Effort_f=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).('Effort_f');   
        
        Effort_e_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Effort_e_Amp');
        Effort_o_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Effort_o_Amp');

        figure(6)
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,Effort_f,'DisplayName',sprintf('E_f, Trial: %d',IndTrials(iTrial)),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        hold on
        plot(FrameInd,Effort_e,'DisplayName',sprintf('E_e, Trial: %d',IndTrials(iTrial)),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,Effort_e+Effort_o,'DisplayName',sprintf('E_e+E_o ,Trial: %d',IndTrials(iTrial)),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
        ylabel('Estimated % Effort')
        xlabel('Frames')


        figure(7)
        subplot(length(TestFolders),1,iTest)
        plot(FrameInd,Effort_f,'DisplayName',sprintf('E_f, Trial: %d',IndTrials(iTrial)),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        hold on
        plot(FrameInd,Effort_e_Amp,'DisplayName',sprintf('E_e, Trial: %d',IndTrials(iTrial)),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,Effort_e_Amp+Effort_o_Amp,'DisplayName',sprintf('E_e+E_o ,Trial: %d',IndTrials(iTrial)),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
        ylabel('Estimated % Effort (Ron s)')
        xlabel('Frames')
    end
    figure(4)
    title(ttl)
    figure(5)
    title(ttl)
    legend('Location','NorthWest')
    grid on
    ylim([-10 100])
end


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
g=strings(1,length(VarNames)-4)+NaN;
x=ones(0,4);

g_dropped=strings(1,length(VarNames)-4)+NaN;
x_dropped=ones(0,4);
x_feat=[];

AnaStruct=sprintf("%s_ana",TestFolders(1));

FeatLabels=string(S.(AnaStruct).AnaPar.FeatLabels);
ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);
sMVCNormLevel=0; % %MVC value for normalization coefficient
clear DrpOrderCat
for iFeat=1:1
    FeatLabel=FeatLabels(iFeat);
    
    for iTest=1:length(TestFolders)
        TestStruct=sprintf("%s_test",TestFolders{iTest});
        AnaStruct=sprintf("%s_ana",TestFolders{iTest});

        RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;

        for iVoli=1:length(VoliMVCLevels)
            sMVC_Vec=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('sMVC');
            vMVC_Vec=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('vMVC');

            Voli_NormCoeff=S.(AnaStruct).(ExpLabel).MAV_Mean_Reps.('Mean')...
                (sMVC_Vec==0 & vMVC_Vec==VoliMVCLevels(iVoli));

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
                    x_dropped=[x_dropped ;x_droppedfeat x_droppednorm x_droppedtarget x_droppedframenum ];

                end
            end
        end
    end
end

Dropped_stats=array2table([[x ;x_dropped] [g(2:end,:); g_dropped(2:end,:)] ],'VariableNames',[VarNames]);

S.(AnaStruct).(ExpLabel).Dropped_stats=Dropped_stats;
writetable( Dropped_stats, 'dropped_stats2.csv')

%% Dropped Frames based Occlusion estimation 
% Occ = 0-stim - avg(dropped frames)

TimeRange=[10 15] ;
AnaLabel=sprintf("%s_ana",TestFolders{1});
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
VoliMVCLevels=[10 20 30 40];
StimMVCLevels=[10 20 30];

Occ=table([],[],[],[],[],[],[],[],[],[],'VariableNames',["Occ_MAV" "NoStim_MAV"...
    "Dropped_MAV" "NormCoef" "Repeat" "Trial" "StimMVC" "VoliMVC" "Test" "NumofDropped"]);
NormMVC=40;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameRange=[TimeRange(1)*stim_freq TimeRange(2)*stim_freq];
    FrameRangeInd=[FrameRange(1): FrameRange(2)];
    NumofDropped=S.(TestLabel).(ExpLabel).num_of_dropped;
    NumofTrial=S.(TestLabel).(ExpLabel).NumofTrials;
    DropOccTest= [[] [] [] [] [] [] [] [] [] [] ];

    for iStim=1:length(StimMVCLevels)
        for iVoli=1:length(VoliMVCLevels)
            
            TrialNums=find_trialnum(VoliMVCLevels(iVoli),StimMVCLevels(iStim),...
                S.(TestLabel).(ExpLabel).RepTableMat);
            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('vMVC');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
            NoStim_mean=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps...
                (vMVC==VoliMVCLevels(iVoli) & sMVC==0,:).('Mean');
            NormCoef=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps...
                (vMVC==NormMVC & sMVC==0,:).('Mean');
            
            for iRep=1:length(TrialNums)
                TrialLabel=sprintf("Trial_%d",TrialNums(iRep));
                
                DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
                MAV_dropped=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat...
                    (FrameRange(1)>=DroppedFrames & DroppedFrames<=FrameRange(2),:).('MAV');
                
                occ(iRep)=(mean(MAV_dropped)-NoStim_mean)/NormCoef;
                NoStim(iRep)=NoStim_mean/NormCoef;
                DroppedMAV(iRep)=mean(MAV_dropped)/NormCoef;
                Repeats(iRep)=sprintf("Rep_%d",iRep);
                Trials(iRep)=TrialLabel;
                StimMVC(iRep)=StimMVCLevels(iStim);
                VoliMVC(iRep)=VoliMVCLevels(iVoli);
                Tests(iRep)=string(TestFolders(iTest));
                NumDropped(iRep)=sprintf("%d_Drops",NumofDropped);
                NormCoefs(iRep)=NormCoef;
            end
            DropOccTest=[DropOccTest; occ' NoStim' DroppedMAV' NormCoefs' Repeats'...
                Trials' StimMVC' VoliMVC' Tests' NumDropped' ];
        end
    end
    DropOccTest=table(DropOccTest(:,1),DropOccTest(:,2),DropOccTest(:,3),...
        DropOccTest(:,4),DropOccTest(:,5),DropOccTest(:,6),DropOccTest(:,7),...
        DropOccTest(:,8),DropOccTest(:,9),DropOccTest(:,10),'VariableNames',...
        ["Occ_MAV" "NoStim_MAV" "Dropped_MAV" "NormCoef" "Repeat" "Trial"...
        "StimMVC" "VoliMVC" "Test" "NumofDropped"]);
    
    S.(AnaLabel).(ExpLabel).DropOccTest=DropOccTest;
    Occ=[Occ; DropOccTest];
end

writetable( Occ, 'dropped_occ.csv')

%% Plotting dropped_occ data 
VoliMVCLevels=[10 20 30 40];
StimMVCLevels=[ 10 20 30];
cm = lines(length(StimMVCLevels)+1);
clc
MeanDropMAV = zeros (length(TestFolders),length(StimMVCLevels),length(VoliMVCLevels));

for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
    
    DropOcc= S.(AnaLabel).(ExpLabel).DropOccTest;
    for iVoli =1:length(VoliMVCLevels)
        VoliVal=sprintf("%d",VoliMVCLevels(iVoli));

        for iStim=1:length(StimMVCLevels)
            StimVal=sprintf("%d",StimMVCLevels(iStim));
            
            ind=DropOcc.('VoliMVC')== VoliVal & DropOcc.('StimMVC')==StimVal;
            MeanDropMAV(iTest,iStim,iVoli)=mean(str2double(DropOcc(ind,:).('Dropped_MAV')));
        end
    end
end

[t,r,c]=size(MeanDropMAV);
clear mC
for iR=1:r
    mC(iR,:)=squeeze(mean(MeanDropMAV(:,iR,:)));
    sdC(iR,:)=squeeze(std(MeanDropMAV(:,iR,:)));

    figure(10)
    plot(VoliMVCLevels,mC(iR,:),'o','LineWidth',1,'Color',cm(iR,:),...
        'DisplayName',sprintf('Stim: %d%%',StimMVCLevels(iR)))
    hold on

    p=polyfit(VoliMVCLevels,mC(iR,:),1);
    plot(VoliMVCLevels,polyval(p,VoliMVCLevels),'LineWidth',2,'Color',...
        cm(iR,:),'DisplayName',sprintf('Stim: %d%%',StimMVCLevels(iR)))
    
end
legend
    







