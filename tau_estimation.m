%% Time constant estimation of RCCurve Trials
% Identify tau (time constant) for a given RCCurve trial
%% Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12"];
TestFiles=["jan7_test","jan11_test" "jan12_test"];
StructstoLoad=["ExpPar","RCCurveTrials"]; 
ExpName=["RCCurveTrials"];
S = load_test(TestFolders,TestFiles,StructstoLoad); % 1- folder name (string), 2- substructures exp numbers

%% Filter Design
fs=S.(TestFiles{1}).ExpPar.fs;
d1 = designfilt("lowpassiir",'FilterOrder',3, ...
    'HalfPowerFrequency',0.01,'DesignMethod',"butter"); %,'SampleRate',fs

% fvtool(d1)
%% Main code 
% This code will execute the estimation of the time constant 
% Redo Trials must be incorporated 

clc
TauTime=.37;
clear Tau F_mat IndTau

for iTest=1:length(TestFiles)
    TestLabel=sprintf("%s",TestFiles{iTest});
    TestStruct=TestLabel;
    NumofTrials=S.(TestLabel).(ExpName).NumofTrials;
    
    DataInd= S.(TestStruct).ExpPar.DataInd;
    DataVars=DataInd.Properties.VariableNames;
    iForce=(S.(TestLabel).ExpPar.DataInd.("Force"));
    iTime=(S.(TestLabel).ExpPar.DataInd.("Time"));
    
    PWVal=S.(TestLabel).(ExpName).PWTrials;
    TurnOffTime=S.(TestLabel).(ExpName).PWProfile(1,4);
    AvgRangePostOff=[TurnOffTime+.5 TurnOffTime+.8];
    AvgRangePreOff=[TurnOffTime-1.5 TurnOffTime];    
    fs=S.(TestFiles{iTest}).ExpPar.fs_calc;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf("Trial_%d",iTrial);

        F=S.(TestLabel).(ExpName).(TrialLabel).data.(DataVars{iForce});
        T=S.(TestLabel).(ExpName).(TrialLabel).data.(DataVars{iTime});
         %% 
         %There is no need for this 
         %%
        F_filtered = filtfilt(d1,F);
        
        Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
        PreAvgForce= mean(F_filtered(Ind));
        Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
        PostAvgForce=mean(F_filtered(Ind));
%         t_end=10.1-TurnOffTime;

        tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1)),-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce);
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
    
    TauVars={'PW','Tau Ind','Time Const' ,'TestNum','TrialNum','PreAvgForce','PostAvgForce'};
    Tau_table=array2table(Taus,'VariableNames',TauVars);
    Redo= false(height(Tau_table),1);
    T_redo = table(Redo,'VariableNames',"Redo");
    Tau_table= [Tau_table T_redo];
    Tau.(TestLabel).(ExpName).TimeCons.TurnOffTime=TurnOffTime;
    Tau.(TestLabel).(ExpName).TimeCons.AvgRangePostOff=AvgRangePostOff;
    Tau.(TestLabel).(ExpName).TimeCons.AvgRangePreOff=AvgRangePreOff;
    Tau.(TestLabel).(ExpName).TimeCons.FiltDesign=d1;
    Tau.(TestLabel).(ExpName).TimeCons.Taus=Taus;
    Tau.(TestLabel).(ExpName).TimeCons.Tau_table=Tau_table;
    
    F_table=array2table(F_mat,'VariableNames',VarNames);
    clear F_mat
    Tau.(TestLabel).(ExpName).TimeCons.F_filtered=F_table;

    %%
    % Similar algo for the redo trials
    %%
    
    clear F_mat VarNames Taus IndTau

    RedoTrials=S.(TestLabel).(ExpName).RedoTrials;
    S.(TestLabel).(ExpName).TimeConsRedo.RedoTrials=RedoTrials;
    PWVal=S.(TestLabel).(ExpName).PWTrials;

    if ~isempty(RedoTrials)
        for iTrial=1:length(RedoTrials)
            TrialLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));

            F=S.(TestLabel).(ExpName).(TrialLabel).data(:,iForce);
            T=S.(TestLabel).(ExpName).(TrialLabel).data(:,iTime);
            
            F_filtered = filtfilt(d1,F);
            Ind= T> AvgRangePreOff(1) & T<AvgRangePreOff(2);
            PreAvgForce= mean(F_filtered(Ind));
            Ind= T> AvgRangePostOff(1) & T<AvgRangePostOff(2);
            PostAvgForce=mean(F_filtered(Ind));
            
%             t_end=mean(AvgRangePreOff)-TurnOffTime;
%             tau=-t_end*log((PreAvgForce-PostAvgForce)/PreAvgForce);
            tau=est_tau(F_filtered(T> TurnOffTime & T<AvgRangePostOff(1)),-TurnOffTime+T(T> TurnOffTime & T<AvgRangePostOff(1)),PreAvgForce);

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
            Taus(iTrial,1)= PWVal(iTrial);

            F_mat(:,iTrial)=F_filtered()';
            VarNames{iTrial}=char(TrialLabel);
        end

        Tau_table_redo=array2table(Taus,'VariableNames',TauVars);
        
        Redo= true(height(Tau_table_redo),1);
        T_redo = table(Redo,'VariableNames',"Redo");
        Tau_table_redo= [Tau_table_redo T_redo];
        Tau.(TestLabel).(ExpName).TimeConsRedo.TurnOffTime=TurnOffTime;
        Tau.(TestLabel).(ExpName).TimeConsRedo.AvgRangePostOff=AvgRangePostOff;
        Tau.(TestLabel).(ExpName).TimeConsRedo.AvgRangePreOff=AvgRangePreOff;
        Tau.(TestLabel).(ExpName).TimeConsRedo.FiltDesign=d1;
        Tau.(TestLabel).(ExpName).TimeConsRedo.Taus=Taus;
        Tau.(TestLabel).(ExpName).TimeConsRedo.Tau_table=Tau_table_redo;

        F_table_redo=array2table(F_mat,'VariableNames',VarNames);
        clear F_mat

        Tau.(TestLabel).(ExpName).TimeConsRedo.F_filtered=F_table_redo;
        F_table=[F_table F_table_redo];
        Tau_table=[Tau_table; Tau_table_redo];        

    end
    Tau.(TestLabel).(ExpName).Tau_table=Tau_table;
    Tau.(TestLabel).(ExpName).F_table=F_table;

end

% Describe with Stats
% Mean, std
 
clear MeanTau StdTau
% Incorporate redos

for iTest= 1:length(TestFiles)
    TestLabel= sprintf("%s_test",TestFolders{iTest});
    
    Tau_table=Tau.(TestLabel).(ExpName).Tau_table;
    RedoInd=table2array(Tau_table(:,"Redo"));
    Redo_table=Tau_table(RedoInd,"TrialNum");
    Tau_table(table2array(Redo_table),:)=Tau_table(table2array(Tau_table(:,"Redo")),:);
    Tau_table(RedoInd,:)=[];
    Tau.(TestLabel).(ExpName).Tau_incorp=Tau_table;
end

Tau_stats=table([],[],[],[],'VariableNames',["PW","mean_Time Const","std_Time Const","Test"]);
% Finding groups
for iTest=1:length(TestFiles)
    TestLabel=sprintf("%s_test",TestFolders{iTest});
    Tau_{iTest}=Tau.(TestLabel).(ExpName).Tau_incorp;
    [ID_PW{iTest},PWPoints{iTest}]=findgroups(Tau_{iTest}(:,1));
%     PWPoints{iTest} = renamevars(PWPoints{iTest},["PW"],[sprintf('%s PW',TestFolders{iTest})]);
    PWPoints{iTest}
end

PWInd{iTest}=ID_PW{iTest} == Ind(iTest);

for iTest=1:length(TestFiles)
    TestLabel=sprintf("%s",TestFolders{iTest});

    for PWPointsInd=1:height(PWPoints{iTest})

        if isempty( PWPointsInd)
           disp('PW value was not found') 
           return
        end

        PWInd{iTest}=ID_PW{iTest} == PWPointsInd;
        MeanTau{iTest}(PWPointsInd,:) = varfun(@mean, Tau_{iTest}(PWInd{iTest},'Time Const'), 'InputVariables', @isnumeric);
        StdTau{iTest}(PWPointsInd,:) = varfun(@std, Tau_{iTest}(PWInd{iTest},'Time Const'), 'InputVariables', @isnumeric);

    end
    Tau_stats=[Tau_stats; PWPoints{iTest} MeanTau{iTest} StdTau{iTest}...
        table([TestLabel+strings(height(PWPoints{iTest}),1)],'VariableNames',"Test")];
end

MTau = array2table(ones(1,length(TauVars)+1),'VariableNames',[TauVars "Redo"]);

for iTest=1:length(TestFolders)
    
    MTau=[MTau; Tau_{iTest}];
end
MTau(1,:)=[]
Tau_stats

% Export Add to Original File

% Export as a mat file, and table

FileNameExtension="_tau";
for iTest=1:length(TestFiles)
    ChrTest=char(TestFiles{iTest});
    TestLabel=sprintf("%s",ChrTest(1:end-5));

    FoldLabel=TestFolders{iTest};
    DirLabel=sprintf('%s/%s%s_test',FoldLabel,FoldLabel,FileNameExtension);
    save(DirLabel,'-struct','Tau',string(ChrTest))
    %%
    %Tau table ".csv" save
    DirLabelCSV=sprintf('%s/%s_test%s.csv',FoldLabel,FoldLabel,FileNameExtension);
    writetable( Tau.(ChrTest).(ExpName).Tau_incorp, DirLabelCSV)
    %%
    %Tau stats ".csv" save
    DirLabelCSV=sprintf('%s/%s_taustats.csv',FoldLabel,FoldLabel);
    writetable( Tau_stats(Tau_stats.Test==FoldLabel,:), DirLabelCSV)
    %all the results in one file
    writetable( MTau,'tau_estimates.csv')

    
end

%% ----------------LoaD-------------------
clear all
% TestFiles=["dec5_tau_test","nov28_2_tau_test","nov27_tau_test","nov8_tau_test"];
% TestFolder=["dec5","nov28_2","nov27","nov8"];
TestFiles=["jan11_tau_test"];
TestFolder=["jan11"];
ExpName=["RCCurveTrials"];
K = load_test(TestFolder,TestFiles);

for iTest=1:length(TestFiles)
    TestLabel=sprintf("%s_test",TestFolder{iTest});
    Tau{iTest}=K.(TestLabel).(ExpName).Tau_table;
    [ID_PW{iTest},PWPoints{iTest}]=findgroups(Tau{iTest}(:,1));
    PWPoints{iTest} = renamevars(PWPoints{iTest},["PW"],[sprintf('%s PW',TestFolder{iTest})]);
    PWPoints{iTest}
    
    %load tables 
    DirLabelCSV=sprintf('%s/%s_taustats.csv',TestFolder{iTest},TestFolder{iTest});
    Tau_table{iTest}= readtable(DirLabelCSV);

end
%% Plotting
PWPoints{1:length(TestFolder)}

close all
PlotPW=[ 120];
for iTest=1:length(TestFolder)
    TestLabel=sprintf('%s_test',TestFolder{iTest});
    
    Ind(iTest)=find( table2array(PWPoints{iTest})==PlotPW(iTest));

    if isempty( Ind(iTest))
       disp('PW value was not found') 
       return
    end

    PWInd{iTest}=ID_PW{iTest} == Ind(iTest);
    Tau{iTest}(PWInd{iTest},:);
    TrialNums{iTest}=Tau{iTest}(:,'TrialNum');

    PWInd_F=logical([1 PWInd{iTest}']');
    F_PW=K.(TestLabel).(ExpName).F_table(:,PWInd_F);
    
    %%
    %Plotting here
    h= figure(1);
    set(h, 'Visible', 'on');
    subplot(length(TestFolder),1,iTest)
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
    str_test=sprintf("Test: %s ",TestFolder{iTest} );
    title(strcat(trials_str1,trials_str2,PW_str,",",str_test))
    lgd_str{1}= sprintf ("Trial %d",TrialVals(1));
    lgd_str{2}= sprintf ("Trial %d",TrialVals(2));
    lgd_str{3}= sprintf ("Trial %d",TrialVals(3));
    lgd_str{4}= sprintf ("Norm PW");

    legend(lgd_str,'Location','NorthWest','AutoUpdate','off')
    xlim([5 11.9])
    grid on
    

    %% 
    %annotations and mark important points
    %%
    a = get(gca,'Children');
    y1data = get(a, 'YData');
    y1min=min( [min(y1data{1}) min(y1data{2}) min(y1data{3}) ]);
    y1max=max( [max(y1data{1}) max(y1data{2}) max(y1data{3}) ]);
    TauInd(:,iTest)=table2array(Tau{iTest}(PWInd{iTest},'Tau Ind'));
    FiltForceVal(:,iTest)=table2array(Tau{iTest}(PWInd{iTest},'FiltForceVal'));
    
    PreAvgForce(:,iTest)=table2array(Tau{iTest}(PWInd{iTest},'PreAvgForce'));
    PostAvgForce(:,iTest)=table2array(Tau{iTest}(PWInd{iTest},'PostAvgForce'));
    AvgRangePreOff=K.(TestLabel).(ExpName).TimeCons.AvgRangePreOff;
    AvgRangePostOff=K.(TestLabel).(ExpName).TimeCons.AvgRangePostOff;
    TurnOffTime=K.(TestLabel).(ExpName).TimeCons.TurnOffTime;
    
    plot([table2array(F_PW(TauInd(:,iTest),1)) table2array(F_PW(TauInd(:,iTest),1))]',...
        [0 FiltForceVal(1,iTest);0 FiltForceVal(2,iTest); 0 FiltForceVal(3,iTest)]','--','Color','k')
    plot([TurnOffTime TurnOffTime]',[0 y1max]','--','Color','k')
    plot([AvgRangePreOff(1) AvgRangePreOff(2)],[PreAvgForce(:,iTest) PreAvgForce(:,iTest)],'--','Color','k')
    plot([AvgRangePostOff(1) AvgRangePostOff(2)],[PostAvgForce(:,iTest) PostAvgForce(:,iTest)],'--','Color','k')
end


