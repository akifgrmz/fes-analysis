%% Occlusion Analysis
% This file contains all the occlusion related analyses
%% Data Inject 
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "apr20" "may19" "oct11" "oct18" "oct25"];
% TestFolders=["jan7" "jan11" "jan12" "apr20" ];
TestFolders=["feb29_24"];



for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%% Plotting the Occ Trials

close all
clc
lbl='Occ';
PlotVoli=3;
PlotStim=2;
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 12 15 ];
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
    MTau.('Test')= renamecats(MTau.('Test'),TestFolders);
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
%%all the results in one file
writetable( MTau,'tau_estimates2.csv')

%%
% # ---------------LoaD-------------------
clc
clear all
% TestFiles=["dec5_tau_test","nov28_2_tau_test","nov27_tau_test","nov8_tau_test"];
% TestFolder=["dec5","nov28_2","nov27","nov8"];

TestFolders=["jan7" "jan11" "jan12"];
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct11" "oct18" "oct25"];
% TestFolders=["jan7" "jan11" "jan12" "apr20" "may19" "oct11" "oct18" "oct25"];

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
    PWPoints{iTest};
    
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

%% Plotting Dropped Frames with Voli only trials
TestFolders=["jan7" "jan11" "jan12"];

cm = lines(length(TestFolders));
sMVC=0;
DroppedsMVC=[10 20 30];
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
    
    sMVCzero_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
    MAVMean_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('Mean');
    vMVCVal_Reps=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('vMVC');
    
    sMVCInd=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC'); 
    MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean.('Mean');

    vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');
    
    Amp_Modul_Mean=S.(AnaLabel).(ExpLabel).Amp_Modul_Mean.('Mean');
    
        % Voli Only trials

    Ind=(sMVCInd==sMVC);
    Ind_Reps=(sMVCzero_Reps==sMVC);
    p1=polyfit(vMVCVal(Ind),MAVMean(Ind),1);
    p2=polyfit(vMVCVal(Ind),Amp_Modul_Mean(Ind),1);
        % Dropped frames 
        
    for iStim=1:length(DroppedsMVC)
        
        DroppedMean=S.(AnaLabel).(ExpLabel).MAV_Mean.('MeanDropped');
        IndDropped=(sMVCInd==DroppedsMVC(iStim));
        p3=polyfit(vMVCVal(IndDropped),DroppedMean(IndDropped),1);
        p4=polyfit(vMVCVal(Ind),Amp_Modul_Mean(Ind),1);

        figure(1)
        subplot(length(DroppedsMVC),1,iStim)
        plot(vMVCVal(Ind),MAVMean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
        hold on
        plot(vMVCVal(Ind),polyval(p1,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))

        plot(vMVCVal(Ind),DroppedMean(IndDropped),'*','LineWidth',1,'Color',cm(iTest,:),'DisplayName',sprintf("%s Dropped",TestFolders(iTest)))
        plot(vMVCVal(Ind),polyval(p3,vMVCVal(Ind)),'--','LineWidth',2,'Color',cm(iTest,:),'DisplayName',sprintf("%s Dropped",TestFolders(iTest)))

        xlabel('% Voli. MVC')
        ylabel('Mean MAV')
        xlim([0 50])
        legend


    %     subplot(2,1,2)
    %     plot(vMVCVal(Ind),Amp_Modul_Mean(Ind),'o','LineWidth',1,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    %     legend
    %     hold on
    %     plot(vMVCVal(Ind),polyval(p2,vMVCVal(Ind)),'LineWidth',2,'Color',cm(iTest,:),'DisplayName',TestFolders(iTest))
    %     ylabel('% Amp modul')
    %     xlabel('% Voli. MVC')
    %     xlim([0 50])

    end
end



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
Occ=[ [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []];

for iTau=1:length(OccRefs)
    OccRef=OccRefs(iTau);

    for iTest=1:length(TestFolders)
        TestLabel=sprintf("%s_test",TestFolders(iTest));
        AnaLabel=sprintf("%s_ana",TestFolders(iTest));
        FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);
        
        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        OccRefMargin= 1/stim_freq ; % in secs
        PreOffMargin=10/stim_freq; % Avg period before stim turning off 
    %-------------- This section should go to main_analysis 
    % MAV_MAX = avg MAV at the last 2 seconds
        AvgTime=2;
        AvgInd=stim_freq*AvgTime-5;
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('MVC');
        FiltLabel="Unfilt";
        MVC=S.(TestLabel).(ExpLabel).MVC;
        for iTrial=1:S.(TestLabel).(ExpLabel).NumofTrials
            TrialLabel=sprintf('Trial_%d',iTrial );
            
            TrialsMAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                .Feats(end-AvgInd:end,:).('MAV_vEMG'));
            TrialsAmp_MAV(iTrial) = mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                .AmpModulFeats(end-AvgInd:end,:).('Amp_MAV_vEMG'));
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

        MAV_max_theo=polyval(poly1,120);
        AmpModul_max_theo=polyval(poly2,120);

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
                
    %---------------
        ExpLabel= S.(AnaLabel).AnaPar.ExpTable.('Occ');
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        TurnOffTime=S.(TestLabel).(ExpLabel).StimProfile;
        RepMatTable=array2table(S.(TestLabel).(ExpLabel).RepTableMat,...
            'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
            "MVC_Voli" "MVC_Stim" "PW","Done"]);
%         S.(TestLabel).(ExpLabel).RepTableMat=RepMatTable;          
        
        
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
            MVC_Voli=RepMatTable(iTrial,:).('MVC_Voli');
            MVC_Stim=RepMatTable(iTrial,:).('MVC_Stim');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('sMVC');
            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps.('vMVC');
            
            MAV_v=S.(AnaLabel).(ExpLabel).MAV_Mean_Reps(MVC_Voli==vMVC & sMVC==0,:).('Mean');
            AmpModul_MAV_v=S.(AnaLabel).(ExpLabel).Amp_Modul_Mean_Reps(MVC_Voli==vMVC & sMVC==0,:).('Mean');

            MAV_max=MAV_max_theo;            %% MAV_MAX replacement 
            AmpModul_MAV_max=AmpModul_max_theo;
            
            DroppedMAV=S.(AnaLabel).(ExpLabel).MAV_Mean(iTrial,:).('MeanDropped');
            vMVC_iTrial=S.(AnaLabel).(ExpLabel).MAV_Mean(iTrial,:).('vMVC');
            vMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');
            sMVCTrials=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC');
            Mean_NoStimMAV=mean(S.(AnaLabel).(ExpLabel).MAV_Mean(vMVCTrials==vMVC_iTrial & sMVCTrials==0,:).('Mean'));
            
            % 1- E_o = E_v-E_d       (Dropped) 
            % 2- E_o = (E_v-E_d)/mvc*100       (Dropped Effort) 
            % 3- E_o = (E_vprime-E_d)/mvc*100   (Hybrid Effort)

            MAV_Dropped= Mean_NoStimMAV-DroppedMAV;
            MAV_Dropped_mvc= MAV_Dropped/MAV_max*100;
            MAV_Hybrid_mvc= F_vprime_mvc-DroppedMAV/MAV_max*100;
            
            Occ=[ Occ;  F_occ  F_vprime F_sprime F_v MAV_Dropped F_target F_occ_mvc F_vprime_mvc...
                F_sprime_mvc F_v_mvc MAV_Dropped_mvc MAV_Hybrid_mvc F_target_mvc TestFolders(iTest) TauLabels(iTau)...
                "Force" "Unfilt" table2array(RepMatTable(iTrial,:) ) iTrial];
            
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
                
                Occ=[ Occ;  MAV_occ  MAV_vprime MAV_sprime MAV_v MAV_Dropped MAV_target MAV_occ_mvc ...
                    MAV_vprime_mvc MAV_sprime_mvc MAV_v_mvc MAV_Dropped_mvc MAV_Hybrid_mvc MAV_target_mvc...
                    TestFolders(iTest) TauLabels(iTau) "MAV" FiltLabel table2array(RepMatTable(iTrial,:)) iTrial];
            
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
                    MAV_Dropped AmpModul_MAV_target AmpModul_MAV_occ_mvc AmpModul_MAV_vprime_mvc...
                    AmpModul_MAV_sprime_mvc AmpModul_MAV_v_mvc MAV_Dropped_mvc MAV_Dropped_mvc ...
                    AmpModul_MAV_target_mvc TestFolders(iTest) TauLabels(iTau) "Amp_Modul" FiltLabel...
                    table2array(RepMatTable(iTrial,:)) iTrial];
            end
        end

        OccTest=Occ(Occ(:,14)==TestFolders(iTest),:); % 14 is where Test variable is located
        OccTest=table(OccTest(:,1),OccTest(:,2),OccTest(:,3),OccTest(:,4),OccTest(:,5), OccTest(:,6), OccTest(:,7),...
            OccTest(:,8),OccTest(:,9),OccTest(:,10),OccTest(:,11),OccTest(:,12),OccTest(:,13),OccTest(:,14),...
            OccTest(:,15),OccTest(:,16),OccTest(:,17),OccTest(:,18),OccTest(:,19),OccTest(:,20),OccTest(:,21),...
            OccTest(:,22),OccTest(:,23),OccTest(:,24),OccTest(:,25),'VariableNames',...
            [ "Occ" "vprime" "sprime" "v" "Occ_Dropped" "Target" "Occ_mvc" "vprime_mvc" "sprime_mvc"...
            "v_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc" "Target_mvc" "Test" "Tau" "Feat" "Filt" "Target_Level"...
            " Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW" "Done" "Trial"]);
        
        
        S.(AnaLabel).(ExpLabel).OccTest=OccTest;
        clear OccTest
    end
end

OccTable=table(Occ(:,1),Occ(:,2),Occ(:,3),Occ(:,4),Occ(:,5), Occ(:,6), Occ(:,7),...
    Occ(:,8),Occ(:,9),Occ(:,10),Occ(:,11),Occ(:,12),Occ(:,13),Occ(:,14),Occ(:,15),...
    Occ(:,16),Occ(:,17),Occ(:,18),Occ(:,19),Occ(:,20),Occ(:,21),Occ(:,22),Occ(:,23),...
    Occ(:,24),Occ(:,25),'VariableNames',...
    [ "Occ" "vprime" "sprime" "v" "Occ_Dropped" "Target" "Occ_mvc" "vprime_mvc" "sprime_mvc"...
    "v_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc" "Target_mvc" "Test" "Tau" "Feat" "Filt" ...
    "Target_Level" "Stim_Force" "Voli_Force" "MVC_Voli" "MVC_Stim" "PW" "Done" "Trial"]);

% writetable( OccTable,'occlusion_v4.csv')

%% Effort Simulation 
%linear modeling for individual occ predictions 
    % Individual occlusion estimation
clear lm_table
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20"];
lm_table=table();
OccType=["Occ_mvc" "Occ_Dropped_mvc" "Occ_Hybrid_mvc"];
FiltTypes=["Unfilt" "GS" "Comb"];
CoefMat=[];
CoefMat2=[];
RowInd=OccTable.('Feat')=="Force" &...
OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau" &...
OccTable.('MVC_Stim')~="0";

for iType=1:length(OccType)
    
    % 1- Linear fitting for Individualized results
    lm_table.Occ=str2double(OccTable(RowInd,:).(OccType(iType)));
    lm_table.LogOcc=abs(log(str2double(OccTable(RowInd,:).(OccType(iType)))));
    lm_table.VoliMVC=str2double(OccTable(RowInd,:).('MVC_Voli'));
    lm_table.StimMVC=str2double(OccTable(RowInd,:).('MVC_Stim'));
    lm_table.Test=categorical((OccTable(RowInd,:).('Test')));
    
    lm_table.Test=reordercats(lm_table.Test,TestFolders);
    mdl = fitlm(lm_table,'Occ~VoliMVC+StimMVC+Test');
    CoefNum=table2array(mdl.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
    CoefCat=table2array(mdl.Coefficients([1,4:end],'Estimate'));
    
    for iTest=1:length(TestFolders)
        TestInd=lm_table.Test==TestFolders(iTest);

        if iTest==1
            Coefs(iTest,:)=[CoefCat(1) CoefNum'];
        else
            Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum'];
        end
        Effort_o=[];
        Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC(TestInd)+Coefs(iTest,3)*lm_table.StimMVC(TestInd))];

        CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) TestFolders(iTest) OccType(iType) boolean(1) boolean(0)],...
            length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd)...
            lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
 
    end
    

    CoefMat=[ CoefMat; Coefs TestFolders' strings(length(TestFolders),1)+OccType(iType) boolean(ones(length(TestFolders),1))  boolean(zeros(length(TestFolders),1))];
    
    mdlLog = fitlm(lm_table,'LogOcc~VoliMVC+StimMVC+Test');
    CoefNum=table2array(mdlLog.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
    CoefCat=table2array(mdlLog.Coefficients([1,4:end],'Estimate'));
%     Effort_o=[];
    for iTest=1:length(TestFolders)
        TestInd=lm_table.Test==TestFolders(iTest);

        if iTest==1
            Coefs(iTest,:)=[CoefCat(1) CoefNum'];
        else
            Coefs(iTest,:)=[CoefCat(1)+CoefCat(iTest) CoefNum'];
        end
        
        Effort_o=[];
        Effort_o=[ Coefs(iTest,1)+(Coefs(iTest,2)*lm_table.VoliMVC(TestInd)+Coefs(iTest,3)*lm_table.StimMVC(TestInd))];
        
        CoefMat2=[ CoefMat2; repmat([Coefs(iTest,:) TestFolders(iTest) OccType(iType) boolean(1) boolean(1)],...
            length(lm_table.Occ(TestInd)),1) Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd)...
            lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
    end
    
    CoefMat=[ CoefMat; Coefs  TestFolders' strings(length(TestFolders),1)+OccType(iType) boolean(ones(length(TestFolders),1)) boolean(ones(length(TestFolders),1))];
    
%     CoefMat2=[ CoefMat2; repmat([Coefs  TestFolders' strings(length(TestFolders),1)+OccType(iType)...
%         boolean(ones(length(TestFolders),1))  boolean(ones(length(TestFolders),1))],length(lm_table.Occ),1)...
%         lm_table.Occ lm_table.LogOcc lm_table.VoliMVC lm_table.StimMVC [1:length(lm_table.Occ)]' ];

    
    
    % 2- Linear fitting for Averaged results
    lm_table_gen=lm_table(lm_table.('Test') =="jan7" | lm_table.('Test') =="jan11" | lm_table.('Test') =="jan12",:);
    mdl_gen = fitlm(lm_table_gen,'Occ~VoliMVC+StimMVC');
    CoefNum_gen=table2array(mdl_gen.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
    CoefCat_gen=table2array(mdl_gen.Coefficients([1,4:end],'Estimate'));
    
%     Effort_o=[];
%     Effort_o=[ Effort_o; CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC+CoefNum_gen(2)*lm_table.StimMVC];
    CoefMat=[ CoefMat;  CoefCat_gen CoefNum_gen' "Avg" OccType(iType) boolean(0) boolean(0)];

    for iTest=1:length(TestFolders)  
        TestInd=lm_table.Test==TestFolders(iTest);
        Effort_o=[];
        Effort_o=CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC(TestInd)+CoefNum_gen(2)*lm_table.StimMVC(TestInd);

        CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' TestFolders(iTest) OccType(iType) boolean(0) boolean(0)],length(lm_table.Occ(TestInd)),1)...
            Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
    end
    
    mdl_genLog = fitlm(lm_table_gen,'LogOcc~VoliMVC+StimMVC');
    CoefNum_gen=table2array(mdl_genLog.Coefficients({'VoliMVC', 'StimMVC'},'Estimate'));
    CoefCat_gen=table2array(mdl_genLog.Coefficients([1,4:end],'Estimate'));
%     
%     Effort_o=[];
%     Effort_o=[ Effort_o; CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC+CoefNum_gen(2)*lm_table.StimMVC];
    
    CoefMat=[ CoefMat;  CoefCat_gen CoefNum_gen' "Avg" OccType(iType) boolean(0) boolean(1)];
    
    for iTest=1:length(TestFolders)  
        TestInd=lm_table.Test==TestFolders(iTest);
        Effort_o=[];
        Effort_o=CoefCat_gen+CoefNum_gen(1)*lm_table.VoliMVC(TestInd)+CoefNum_gen(2)*lm_table.StimMVC(TestInd);

        CoefMat2=[ CoefMat2; repmat([CoefCat_gen CoefNum_gen' TestFolders(iTest) OccType(iType) boolean(0) boolean(1)],length(lm_table.Occ(TestInd)),1)...
            Effort_o lm_table.Occ(TestInd) lm_table.LogOcc(TestInd) lm_table.VoliMVC(TestInd) lm_table.StimMVC(TestInd) [1:length(lm_table.Occ(TestInd))]' ];
    end
    
end

OccCoef=table(CoefMat(:,1),CoefMat(:,2),CoefMat(:,3),CoefMat(:,4),CoefMat(:,5),CoefMat(:,6),CoefMat(:,7),...
    'VariableNames',["Coeff1" "Coeff2" "Coeff3" "Test" "Type" "Indiv" "Log"]);

OccCoef2=table(CoefMat2(:,1),CoefMat2(:,2),CoefMat2(:,3),CoefMat2(:,4),CoefMat2(:,5),CoefMat2(:,6),CoefMat2(:,7),...
    CoefMat2(:,8),CoefMat2(:,9),CoefMat2(:,10),CoefMat2(:,11),CoefMat2(:,12),CoefMat2(:,13),...
    'VariableNames',["Coeff1" "Coeff2" "Coeff3" "Test" "Type" "Indiv" "Log" "Effort_o" "Occ" "LogOcc" "MVC_Voli" "MVC_Stim" "Trial"]);

% writetable( OccCoef2,'occ_coef.csv')
%% Plotting fittings 



%%
    %Average Occlusion estimation
% OccCoef=[-1.7 -0.109 0.919];  % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fs_primeCoef=[2.117 0.061 0.118 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]
% Fv_primeCoef=[-1.361 0.891 0.848 ]; % These coefficienst are from Rstudio results: [resid Voli Stim]
%
% Implementing the correction for occlusion
ExpstoAna= ["Occ"];
AnaLabel=sprintf('%s_ana',TestFolders(1));
OccKickOffLevel=0.3;  % Percent effort
for iExp=1:length(ExpstoAna)
    ExpLabels(iExp)=S.(AnaLabel).AnaPar.ExpTable.(ExpstoAna(iExp));
end

ErrMat=[ [] [] [] [] [] [] [] [] [] [] [] [] []];
LogModel=["true" "false"];
FeattoAna=["Filt_MAV_vEMG"];
RampTime=[5 10];
ConstTime=[10 15];
TurnOffTime=[15 17];
RampInd=[stim_freq*RampTime(1):stim_freq*RampTime(2)];
ConstInd=[stim_freq*ConstTime(1):stim_freq*ConstTime(2)];
FrameInd=[ RampInd ConstInd];
clear Effort
for iModel=1:length(LogModel)
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
    
    RepTableMat=array2table(S.(TestLabel).(ExpLabel).RepTableMat,...
        'VariableNames',["Target_Level" "Stim_Force" "Voli_Force"...
        "MVC_Voli" "MVC_Stim" "PW","Done"]);
    
    for iExp=1:length(ExpLabels)
        ExpLabel=ExpLabels(iExp);
        
        NumofTrials=S.(TestLabel).(ExpLabel).NumofTrials;
        for iTrial=1:NumofTrials 
            TrialLabel=sprintf("Trial_%d",iTrial);
            
            DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
            PWofFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).PWofFrames;
            StimMVC=S.(TestLabel).(ExpLabel).RepTableMat(iTrial,5);
            StimMVCofFrames=gompertz(RCVar,PWofFrames)/MVC*100+0.00001;

            KeepInd=S.(AnaLabel).(ExpLabel).(TrialLabel).KeepInd;

            % 1- Force 
            Force=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).ForceFrames(:,KeepInd))';
            Effort_f=Force/F_MAX*100;
            TargetFrames=mean(S.(AnaLabel).(ExpLabel).TargetFrames(:,KeepInd))';
            Target_mvc=RepTableMat(iTrial,:).('Target_Level')/F_MAX*TargetFrames;
            VoliEffort=RepTableMat(iTrial,:).('MVC_Voli')*TargetFrames;
            StimEffort=RepTableMat(iTrial,:).('MVC_Stim')*TargetFrames;
            
            IndTrials=find_trialnum(RepTableMat(iTrial,:).('MVC_Voli'), ...
                RepTableMat(iTrial,:).('MVC_Stim'),S.(TestLabel).(ExpLabel).RepTableMat);
            
%             for iFrame=1:length(Effort_f)
%                 Effort_fs_prime(iFrame,1)=Fs_primeCoef(1)+Fs_primeCoef(2)*VoliEffort(iFrame)+Fs_primeCoef(3)*StimEffort(iFrame);
%                 Effort_fv_prime(iFrame,1)=Fv_primeCoef(1)+Fv_primeCoef(2)*VoliEffort(iFrame)+Fv_primeCoef(3)*StimEffort(iFrame);
%             end

            RowIndTemp=OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')=="Force" & ...
            OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau" & OccTable.('Trial')==string(IndTrials');
        
            [~,c]=size(RowIndTemp);
            RowInd=RowIndTemp(:,1);
            for i=1:c-1
                RowInd=RowInd | RowIndTemp(:,i+1);
            end
            %sum(RowInd)

            Fv_prime=str2double(OccTable(RowInd,:).('vprime_mvc'));
%                 Effort_Fv_prime=[ linspace(0,0,stim_freq*5)' ; linspace(0, mean(Fv_prime),stim_freq*5)'; linspace(mean(Fv_prime), mean(Fv_prime),stim_freq*7+1)' ] ;
            Effort_Fv_prime=TargetFrames*mean(Fv_prime);
                
            for iType=1:length(OccType)
                TypeLabel=OccType(iType);
                
                RowInd=OccCoef.('Type') == OccType(iType) & OccCoef.('Test') == TestFolders(iTest) & OccCoef.('Log') == LogModel(iModel);
                CoefIndv=str2double([ OccCoef(RowInd,:).('Coeff1') OccCoef(RowInd,:).('Coeff2') OccCoef(RowInd,:).('Coeff3') ]);

                RowInd=OccCoef.('Type') == OccType(iType) & OccCoef.('Test') == "Avg" & OccCoef.('Log') == LogModel(iModel);
                CoefAvg=str2double([ OccCoef(RowInd,:).('Coeff1') OccCoef(RowInd,:).('Coeff2') OccCoef(RowInd,:).('Coeff3') ]);

                Effort_o=CoefAvg(1)+CoefAvg(2).*VoliEffort+CoefAvg(3).*StimEffort;
                Effort_o_Indv=CoefIndv(1)+CoefIndv(2).*VoliEffort+CoefIndv(3).*StimEffort;
                
                if LogModel(iModel)=="true"
                    Effort_o=exp(Effort_o);
                    Effort_o_Indv=exp(Effort_o_Indv);
                end
                
                for iFilt=1:length(FiltLabels) 
                    FiltLabel=FiltLabels(iFilt);

                    Filt_Feats=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats.('Filt_MAV_vEMG');
                    MAV_Noise=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).MAV_Noise;
                        % Need error against dropped frames 
                    % 2-MAV
%                     Effort_e_MAV=(Filt_Feats-MAV_Noise)/(MAV_MAX-MAV_Noise)*100;
                    Effort_e_MAV=Filt_Feats/MAV_MAX*100;
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
                    
                    Effort_e_Amp=(AmpModul_Feats-AmpModul_MAV_Noise)/(AmpModul_max_theo-AmpModul_MAV_Noise)*100;
                    Effort_Amp=Effort_e_Amp+Effort_o;
                    Effort_Amp_Indv=Effort_e_Amp+Effort_o_Indv;
                    
                    Effort_e_Amp_err=Effort_e_Amp-Effort_Fv_prime;
                    Effort_Amp_err=Effort_Amp-Effort_Fv_prime;
                    Effort_Amp_Indv_err=Effort_Amp_Indv-Effort_Fv_prime;
                    

                    ErrMat=[ErrMat; mean(Effort_e_MAV(ConstInd)) mean(Effort_MAV(ConstInd)) mean(Effort_MAV_Indv(ConstInd))...
                        mean(Effort_e_Amp(ConstInd)) mean(Effort_Amp(ConstInd)) mean(Effort_Amp_Indv(ConstInd))...
                        mean(Effort_e_MAV_err(ConstInd)) mean(Effort_MAV_err(ConstInd)) mean(Effort_MAV_Indv_err(ConstInd))...
                        mean(Effort_e_Amp_err(ConstInd)) mean(Effort_Amp_err(ConstInd)) mean(Effort_Amp_Indv_err(ConstInd))...
                        TypeLabel FiltLabel TrialLabel ExpLabel TestFolders(iTest) LogModel(iModel)...
                        S.(TestLabel).(ExpLabel).RepTableMat(iTrial,:) ];
% 
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
end
%
Effort_Est=table(ErrMat(:,1),ErrMat(:,2),ErrMat(:,3),ErrMat(:,4),ErrMat(:,5), ErrMat(:,6), ErrMat(:,7),...
    ErrMat(:,8),ErrMat(:,9),ErrMat(:,10),ErrMat(:,11),ErrMat(:,12),ErrMat(:,13),ErrMat(:,14),ErrMat(:,15),...
    ErrMat(:,16),ErrMat(:,17),ErrMat(:,18),ErrMat(:,19),ErrMat(:,20),ErrMat(:,21),ErrMat(:,22),ErrMat(:,23),...
    ErrMat(:,24),ErrMat(:,25),'VariableNames',...
    ["Effort_e_MAV" "Effort_MAV" "Effort_MAV_Indv"...
    "Effort_e_Amp" "Effort_Amp" "Effort_Amp_Indv"...
    "Effort_e_MAV_err" "Effort_MAV_err" "Effort_MAV_Indv_err"...
    "Effort_e_Amp_err" "Effort_Amp_err" "Effort_Amp_Indv_err"...
    "Occ_Type" "Filt" "Trial" "Exp" "Test" "LogModel" "Target_Level" "Stim_Force"...
    "Voli_Force" "VoliMVC" "StimMVC" "PW" "Done"]);

% writetable( Effort_Est,'occ_est_error.csv')


%% Plotting Effort
close all
clc
lbl='Occ';
PlotVoli=4;
PlotStim=4;
% TestFolders=["jan7" "jan11" "jan12"];
TimeRange=[1 15];
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 ];
sMVC=StimMVCLevels(PlotStim);
vMVC=VoliMVCLevels(PlotVoli);
cm=lines(6);
ln=["--" ":" "-."];
FiltLabel="GS";

for iType=1:length(OccType)
    TypeLabel=OccType(iType);

    for iTest=1:length(TestFolders)
        AnaLabel=sprintf('%s_ana',TestFolders(iTest));
        TestLabel=sprintf('%s_test',TestFolders(iTest));
        stim_freq=S.(TestLabel).ExpPar.stim_freq;
        FrameInd=[stim_freq*TimeRange(1): stim_freq*TimeRange(2)];

        ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(lbl);
        RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;

        IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
        ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);

        RowInd=OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')=="Force" &...
            OccTable.('MVC_Voli')==num2str(vMVC) & OccTable.('MVC_Stim')==num2str(sMVC) &...
            OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau";
        Fv_prime=str2double(OccTable(RowInd,:).('vprime_mvc'));

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

            % With individual maintained levels
    %         figure(3)
    %         subplot(length(TestFolders),1,iTest)
    %         plot(FrameInd,Effort_f,'DisplayName',sprintf('E_f, Trial: %d',IndTrials(iTrial)),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
    %         hold on
    %         plot(FrameInd,Effort_e,'DisplayName',sprintf('E_e, Trial: %d',IndTrials(iTrial)),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
    %         plot(FrameInd,Effort_e+Effort_o,'DisplayName',sprintf('E_e+E_o ,Trial: %d',IndTrials(iTrial)),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
    %         plot([10 15]*stim_freq, [mean(Fv_prime) mean(Fv_prime)],'DisplayName',sprintf('Fv_{prime}'),'Color',cm(4,:))
    %         ylabel('Estimated % Effort')
    %         xlabel('Frames')
    %         grid on

    % %         figure(4)
    % %         subplot(length(TestFolders),1,iTest)
    % %         plot(FrameInd,Effort_f,'DisplayName',sprintf('E_f'),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
    % %         hold on
    % %         plot(FrameInd,Effort_e_Amp,'DisplayName',sprintf('E_{MAV}'),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
    % %         plot(FrameInd,(Effort_e_Amp+Effort_o_Amp_Indv),'DisplayName',sprintf('E_{MAV}+E_o'),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
    % %         plot([10 15]*stim_freq, [mean(Fv_prime) mean(Fv_prime)],'DisplayName',sprintf('Fv_{prime}'),'Color',cm(4,:))
    % % 
    % %         ylabel('Estimated % Effort')
    % %         xlabel('Frames')
    % %         grid on


    %     figure(1)
    %     title(ttl)


    % %     figure(2)
    % %     title(ttl)
    % %     ylim([-10 80])
    % %         figure(4)
    % %     title(ttl)
    % %     ylim([-10 80])

    end
    figure(1)
    legend('Location','NorthWest')

    
end
%% Plotting for Conf Paper:EMBC
close all
clc


iRep=3;
TestFolders=["jan7"];
TimeRange=[1 15];
sMVC=30; % [ 0 10 20 30 ] 
vMVC=40; % [10 20 30 40 ]
% VoliMVCLevels=[10 20 30 40 ];
% StimMVCLevels=;

cm=lines(6);
ln=["--" ":" "-."];
FiltLabel="GS";
lbl='Occ';

for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=[stim_freq*TimeRange(1): stim_freq*TimeRange(2)];

    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(lbl);
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;

    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
    ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);

    RowInd=OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')=="Force" &...
        OccTable.('MVC_Voli')==num2str(vMVC) & OccTable.('MVC_Stim')==num2str(sMVC) &...
        OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau";
    Fv_prime=str2double(OccTable(RowInd,:).('vprime_mvc'));

    TrialLabel=sprintf('Trial_%d',IndTrials(iRep));

    TypeLabel=OccType(1);
    E_Amp_Force=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Amp');
    TypeLabel=OccType(2);
    E_Amp_Dropped=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Amp');

    E_e_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_e_Amp');
    E_f=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_f');
    E_fv_prime=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Fv_prime');  

    figure(1)
%         subplot(length(OccType),length(TestFolders),iTest+(iType*length(TestFolders)-length(TestFolders)))
    plot(FrameInd/stim_freq,E_e_Amp,'DisplayName',sprintf('E_e'),'Color',cm(2,:),'LineWidth',2)%,'LineStyle',ln(iTrial))
    hold on
    plot([10 16],E_fv_prime(end)*[1 1],'DisplayName',sprintf('V^{''}'),'Color',cm(5,:),'LineWidth',2,'LineStyle',ln(1))

    plot(FrameInd/stim_freq,E_Amp_Force,'DisplayName',sprintf('E_e + E_o'),'Color',cm(3,:),'LineWidth',2)%,'LineStyle',ln(iTrial))
    plot(FrameInd/stim_freq,E_f,'DisplayName',sprintf('E_f'),'Color',cm(1,:),'LineWidth',2)%,'LineStyle',ln(iTrial))

    plot([0 175 350 525]/stim_freq,(sMVC+vMVC)*[0 0 1 1],'DisplayName',sprintf('Target Line'),'Color',cm(6,:),'LineWidth',2)%,'LineStyle',ln(iTrial))

    %     plot(FrameInd/stim_freq,E_Amp_Dropped,'DisplayName',sprintf('E_e + E^d_o'),'Color',cm(4,:),'LineWidth',2)%,'LineStyle',ln(iTrial))

%         plot(FrameInd,mean(E_Amp_Indv,2),'DisplayName',sprintf('E_e+E_{o, indv} ,Trial: %d, %d, %d',IndTrials),'Color',cm(5,:))%,'LineStyle',ln(iTrial))
    title(TestFolders(iTest),'FontSize', 12);
    ylabel('Effort (%)','FontSize', 12);
    xlabel('Time (s)','FontSize', 12);
    grid on


end
figure(1)
legend('Location','NorthWest','FontSize', 12);


%% PLotting Effort Estimation Error 

close all
clc
lbl='Occ';
PlotVoli=1;
PlotStim=4;
% TestFolders=["jan7" "jan12" "apr20"];
TimeRange=[1 15];
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 ];
sMVC=StimMVCLevels(PlotStim);
vMVC=VoliMVCLevels(PlotVoli);
cm=lines(5);
ln=["--" ":" "-."];
FiltLabel="GS";
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=[stim_freq*TimeRange(1): stim_freq*TimeRange(2)];
    
    ExpLabel=S.(AnaLabel).AnaPar.ExpTable.(lbl);
    RepTableMat=S.(TestLabel).(ExpLabel).RepTableMat;
    
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);
    ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);

    RowInd=OccTable.('Test')==TestFolders(iTest) & OccTable.('Feat')=="Force" &...
        OccTable.('MVC_Voli')==num2str(vMVC) & OccTable.('MVC_Stim')==num2str(sMVC) &...
        OccTable.('Filt')=="Unfilt" & OccTable.('Tau')=="3*tau";
    Fv_prime=str2double(OccTable(RowInd,:).('vprime_mvc'));

    for iTrial=1:length(IndTrials)
    
        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        Effort_f=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).('Effort_f');   
        Effort_fv_prime=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).FatFeats(FrameInd,:).('Effort_fv_prime');   
        
        Effort_e_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Effort_e_Amp');
        Effort_o_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats(FrameInd,:).('Effort_o_Amp');
        
        Effort_err_Amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('Effort_err_Amp');
        stdError_amp=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).AmpModulFeats.('stdError_amp');

        figure(10+iTest)
        subplot(2,1,1)
        plot(FrameInd,Effort_f,'DisplayName',sprintf('E_f'),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        hold on
        plot(FrameInd,Effort_e_Amp,'DisplayName',sprintf('E_{MAV}'),'Color',cm(2,:))%,'LineStyle',ln(iTrial))
        plot(FrameInd,(Effort_e_Amp+Effort_o_Amp),'DisplayName',sprintf('E_{MAV}+E_o'),'Color',cm(3,:))%,'LineStyle',ln(iTrial))
        plot([10 15]*stim_freq, [mean(Fv_prime) mean(Fv_prime)],'DisplayName',sprintf('Fv_{prime}'),'Color',cm(4,:))
        plot(FrameInd,Effort_err_Amp(FrameInd),'DisplayName',sprintf('E_err (%s)',TrialLabel),'Color',cm(5,:))%,'LineStyle',ln(iTrial))

        ylabel('Estimated % Effort')
        xlabel('Frames')
        grid on
        
        subplot(2,1,2)
        plot(FrameInd,stdError_amp(FrameInd),'DisplayName',sprintf('E_ste (%s)',TrialLabel),'Color',cm(1,:))%,'LineStyle',ln(iTrial))
        ylim([-2 2])
        hold on
        ylabel('Standard Error')
        xlabel('Frames')
        grid on

    end 
end

        
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
    
        clear OccForce OccMAV OccAmpModul

    for iVoli=1:length(VoliLevels)
        VoliLabel=VoliLevels(iVoli);
        
        for iStim=2:length(StimLevels)
            StimLabel=StimLevels(iStim);
            
            % 1- Force
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Force'...
                &  OccTest.('Tau')=='3*tau';
            
            OccForce(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));
            OccForceMap(iTest,iVoli,iStim-1)=OccForce(iVoli,iStim-1);
            
            % 2- MAV
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='MAV'...
                &  OccTest.('Tau')=='3*tau';
            
            OccMAV(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));
            OccMAVMap(iTest,iVoli,iStim-1)=OccMAV(iVoli,iStim-1);
            
            % 3- Amplitude Modul.
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Amp_Modul'...
                &  OccTest.('Tau')=='3*tau';
            
            OccAmpModul(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('Occ_mvc')));
            OccAmpMap(iTest,iVoli,iStim-1)=OccAmpModul(iVoli,iStim-1);
            
            % 4- F_sprime_mvc
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Force'...
                &  OccTest.('Tau')=='3*tau';
            
            Fs_prime(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('sprime_mvc')));
            Fs_primeMap(iTest,iVoli,iStim-1)=Fs_prime(iVoli,iStim-1);
            
            % 5- Mean F_fprime
            RowInd=OccTest.('MVC_Stim')==StimLabel & OccTest.('MVC_Voli')==VoliLabel...
                & OccTest.('Filt')=='Unfilt' &  OccTest.('Feat')=='Force'...
                &  OccTest.('Tau')=='3*tau';
            
            Fv_prime(iVoli,iStim-1)=mean(str2double(OccTest(RowInd,:).('vprime_mvc')));
            Fv_primeMap(iTest,iVoli,iStim-1)=Fs_prime(iVoli,iStim-1);
            
        end
    end
    
    dx=1;
    VoliLevels=str2double(VoliLevels);
    StimLevels=str2double(StimLevels(2:end));
    vx=[VoliLevels(1):1/(StimLevels(end)-StimLevels(1))/dx:VoliLevels(end)];
    vy=[StimLevels(1):1/(VoliLevels(end)-VoliLevels(1))/dx:StimLevels(end)];

    vq_force=interpn(VoliLevels,StimLevels,OccForce,vx',vy);
    vq_mav=interpn(VoliLevels,StimLevels,OccMAV,vx',vy);
    vq_amp=interpn(VoliLevels,StimLevels,OccAmpModul,vx',vy);

    S.(AnaLabel).(ExpLabel).OccForce_Base=OccForce;
    S.(AnaLabel).(ExpLabel).OccForce_Interp=vq_force;

    S.(AnaLabel).(ExpLabel).OccMAV_Base=OccMAV;
    S.(AnaLabel).(ExpLabel).OccMAV_Interp=vq_mav;

    S.(AnaLabel).(ExpLabel).OccAmpModul_Base=OccAmpModul;
    S.(AnaLabel).(ExpLabel).OccAmpModul_Interp=vq_amp;

    S.(AnaLabel).(ExpLabel).Interp_vx=vx;
    S.(AnaLabel).(ExpLabel).Interp_vy=vy;

end

dx=1;
Mean_OccForceMap=squeeze(mean(OccForceMap));
vx=[VoliLevels(1):1/(StimLevels(end)-StimLevels(1))/dx:VoliLevels(end)];
vy=[StimLevels(1):1/(VoliLevels(end)-VoliLevels(1))/dx:StimLevels(end)];

vq=interpn(VoliLevels,StimLevels,Mean_OccForceMap,vx',vy);

% figure
% mesh(vx,vy,vq)
% xlabel('Stim MVC(%)')
% ylabel('Voli MVC(%)')
% zlabel('Occlusion MVC (%)')
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

Dropped_stats=array2table([[x ;x_dropped] [g(2:end,:); g_dropped(2:end,:)] ],'VariableNames',[VarNames]);

S.(AnaStruct).(ExpLabel).Dropped_stats=Dropped_stats;
writetable( Dropped_stats, 'dropped_stats2.csv')

%% Dropped Frames based Occlusion Estimation 
% Occ = 0-stim - avg(dropped frames)
% TestFolders=["jan7" "jan11" "apr20" "mar16"];


TimeRange=[10 15] ;
AnaLabel=sprintf("%s_ana",TestFolders{1});
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
VoliMVCLevels=[10 20 30 40];
StimMVCLevels=[10 20 30];

FiltLabels=string(S.(AnaLabel).AnaPar.FiltLabels);

Occ=table([],[],[],[],[],[],[],[],[],'VariableNames',["MAV" "MAV_Type"...
    "NormCoef" "Repeat" "Trial" "StimMVC" "VoliMVC" "Test" "NumofDropped"]);
NormMVC=40;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf("%s_ana",TestFolders{iTest});
    TestLabel=sprintf("%s_test",TestFolders{iTest});

    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameRange=[TimeRange(1)*stim_freq TimeRange(2)*stim_freq];
    FrameRangeInd=[FrameRange(1): FrameRange(2)];
    NumofDropped=S.(TestLabel).(ExpLabel).num_of_dropped;
    NumofTrial=S.(TestLabel).(ExpLabel).NumofTrials;
    
    sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC'); 
    MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean.('Mean');
    vMVCVal=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');
        
    Ind=(sMVCzero==0);
    Ind_Reps=(sMVCzero_Reps==0);
    p=polyfit(vMVCVal(Ind),MAVMean(Ind),1);
    NormCoef=polyval(p,100);
    
    DropOccTest= [[] [] [] [] [] [] [] []];
    
    for iStim=1:length(StimMVCLevels)
        for iVoli=1:length(VoliMVCLevels)
            
            TrialNums=find_trialnum(VoliMVCLevels(iVoli),StimMVCLevels(iStim),... 
                S.(TestLabel).(ExpLabel).RepTableMat);
            
            vMVC=S.(AnaLabel).(ExpLabel).MAV_Mean.('vMVC');
            sMVC=S.(AnaLabel).(ExpLabel).MAV_Mean.('sMVC');
            NoStim=S.(AnaLabel).(ExpLabel).MAV_Mean(vMVC==VoliMVCLevels(iVoli) & sMVC==0,:).('Mean');
            
            for iRep=1:length(TrialNums)
                TrialLabel=sprintf("Trial_%d",TrialNums(iRep));
                
                
                DroppedFrames=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames;
                MAV_dropped=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat...
                    (FrameRange(1)>=DroppedFrames & DroppedFrames<=FrameRange(2),:).('MAV');
                
                occMAV(iRep)=-(mean(MAV_dropped)-mean(NoStim))/NormCoef;
                NoStimMAV(iRep)=NoStim(iRep)/NormCoef;
                DroppedMAV(iRep)=mean(MAV_dropped)/NormCoef;
                Repeats(iRep)=sprintf("Rep_%d",iRep);
                Trials(iRep)=TrialLabel;
                StimMVC(iRep)=StimMVCLevels(iStim);
                VoliMVC(iRep)=VoliMVCLevels(iVoli);
                Tests(iRep)=string(TestFolders(iTest));
                NumDropped(iRep)=sprintf("%d_Drops",NumofDropped);
                NormCoefs(iRep)=NormCoef;
                
                DropOccTest=[DropOccTest; [occMAV(iRep) NoStimMAV(iRep) DroppedMAV(iRep)]'...
                    ["Occ" "NoStim" "Dropped"]' ones(3,1)*NormCoefs(iRep) strings(3,1)+Repeats(iRep)...
                    strings(3,1)+Trials(iRep) ones(3,1)*StimMVC(iRep) ones(3,1)*VoliMVC(iRep)...
                    strings(3,1)+Tests(iRep) strings(3,1)+NumDropped(iRep) ];
                
                
                for iFilt=1:length(FiltLabels)
                    FiltLabel=FiltLabels(iFilt);

                    MAV_Filt=mean(S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel)...
                        .Feats(FrameRange(1):FrameRange(2),:).('MAV_vEMG'))/NormCoef;

                    DropOccTest=[DropOccTest; MAV_Filt...
                        FiltLabel NormCoefs(iRep) Repeats(iRep)...
                        Trials(iRep) StimMVC(iRep) VoliMVC(iRep)...
                        Tests(iRep) NumDropped(iRep) ];
                end
                
                
                
            end
        end
    end
    
    DropOccTest=table(DropOccTest(:,1),DropOccTest(:,2),DropOccTest(:,3),...
        DropOccTest(:,4),DropOccTest(:,5),DropOccTest(:,6),DropOccTest(:,7),...
        DropOccTest(:,8),DropOccTest(:,9),'VariableNames',...
        ["MAV" "MAV_Type"  "NormCoef" "Repeat" "Trial"...
        "StimMVC" "VoliMVC" "Test" "NumofDropped"]);
    
    S.(AnaLabel).(ExpLabel).DropOccTest=DropOccTest;
    Occ=[Occ; DropOccTest];
end

writetable( Occ, 'dropped_occ.csv')

%% Effort Correction using Dropped Frames Occ Estimations 
% 1- E= E_e+E_vprime-E_d   (Hybrid)
% 2- E= E_e-E_d+E_v       (Dropped) 



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
    







