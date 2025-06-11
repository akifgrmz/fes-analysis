%effortFES



%%
for iTest=1:length(TestFolders)
    TestFolder=TestFolders(iTest);
    S=tidy_data(TestFolder);
    TestFile=sprintf("%s_test",TestFolder);
    str=sprintf('%s/%s',TestFolder,TestFile);
    save(str,'-struct','S',TestFile)
    str=sprintf("%s.mat is created",TestFile);
    disp(str)
end

%% Data Inject
clear all

TestFolders=["nov21_24" "nov27_24" "dec12_24"];
FileNames=strings(1,length(TestFolders))+"expsave";

for iFile=1:length(TestFolders)
    ExpStruct=sprintf('%s',char(FileNames(iFile)));
    str=sprintf('%s/%s',char(TestFolders(iFile)),ExpStruct);
    M=load (str);
    FldNames = fieldnames(M);

    for iField=1:length(FldNames)
        S.(TestFolders(iFile))=M.(FldNames{iField});
    end
end

%% Initial Plotting
close all
Trials=[12];
for iTest=1:length(TestFolders)
    TestLabel=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpLabels(6);
    EffortFESTable=S.(TestLabel).(ExpLabel).EffortFESTable;
    CompTrialNum=sum(EffortFESTable(:,end));
    TrialNum=length(EffortFESTable(:,end));

    for iTrial=1:length(Trials)
        TrialLabel=sprintf("Trial_%d",Trials(iTrial));
        
         % rows are frames, columns are samples in a frame
        ign=9;  
        FiltvEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltvEMG_f(ign:end,:);
        FiltmWaves_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltmWaves_f(ign:end,:);
        Trig_f=S.(TestLabel).(ExpLabel).(TrialLabel).Trig_f(ign:end,:);
        BlankedEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).BlankedEMG_f(ign:end,:);
        Dropped_f=S.(TestLabel).(ExpLabel).(TrialLabel).DroppedEMG_f(ign:end,:);
        DroppedTrig_f=S.(TestLabel).(ExpLabel).(TrialLabel).DrpdTrig_f(ign:end,:);
        PW_f=S.(TestLabel).(ExpLabel).(TrialLabel).PW_f(ign:end,:);
        
        Time_f=S.(TestLabel).(ExpLabel).(TrialLabel).Time_f(ign:end,:);                
        Tm_f=S.(TestLabel).(ExpLabel).(TrialLabel).Tm_f(ign:end,:);
        Time_f=Tm_f;
        [r,c]=size(Time_f);
        
        FiltvEMG_f=reshape(FiltvEMG_f',1,[]);
        Time_f=reshape(Time_f',1,[]);
        FiltmWaves_f=reshape(FiltmWaves_f',1,[]);
        Trig_f=reshape(Trig_f',1,[]);
        BlankedEMG_f=reshape(BlankedEMG_f',1,[]);                
        Dropped_f=reshape(Dropped_f',1,[]);
        DroppedTrig_f=reshape(DroppedTrig_f',1,[]);
        PW_f=reshape(PW_f',1,[]);

        Hands=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,4:5);
        PW=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,6);
        MAVs=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,1+[16:20]);
        Stim=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,22);
        Target=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,23);

        Time=S.(TestLabel).(ExpLabel).(TrialLabel).time;
        
        ind=100:r*c-100;
        f=figure(iTest);
        % f.Position = [800 300 800 1500 ];    
        % subplot(length(Trials),2,2*iTrial)
        subplot(2,1,2)        
        % plot(Time_f(ind),Trig_f(ind)/200    ,'LineWidth',0.5,"DisplayName","Trig") 
        hold on        
        plot(Time_f(ind),BlankedEMG_f(ind),'LineWidth',1,"DisplayName","BlankedEMG")  
        plot(Time_f(ind),FiltmWaves_f(ind),'LineWidth',1,"DisplayName","FiltmWaves")        
        plot(Time_f(ind),FiltvEMG_f(ind),'LineWidth',1,"DisplayName","FiltvEMG")
        % plot(Time_f(ind),DroppedTrig_f(ind)/100,'LineWidth',0.5,"DisplayName","DroppedTrig")
        % plot(Time_f(ind),Dropped_f(ind),'LineWidth',1,"DisplayName","DroppedEMG")
%         plot(Time_f(ind),PW_f(ind),'LineWidth',1,"DisplayName","PW")
%         legend({'Trigger(a.u.)','Raw EMG (a.u.)'})
        legend
        title(TrialLabel)
        xlabel('Time (s)')
        ylim([-.1 .1])
        
        
        % subplot(length(Trials),2,2*iTrial-1)
        subplot(2,1,1)        
        plot(Time,Hands(:,1),'LineWidth',1,"DisplayName","Paretic")
        hold on
        plot(Time,Hands(:,2),'LineWidth',1,"DisplayName","NonParetic")

        plot(Time,PW,'LineWidth',1,"DisplayName","PW")
        plot(Time,MAVs(:,1),'LineWidth',1,"DisplayName","vEMG MAV")
        % plot(Time,MAVs(:,2),'LineWidth',1,"DisplayName","mWave MAV")
        % plot(Time,MAVs(:,3),'LineWidth',1,"DisplayName","Unfiltered")
        % plot(Time,MAVs(:,4),'LineWidth',1,"DisplayName","Dropped")                
        plot(Time,MAVs(:,5),'LineWidth',1,"DisplayName","Target")
        plot(Time,Stim,'LineWidth',1,"DisplayName","StimCont")
        plot(Time,Target,'LineWidth',1,"DisplayName","Adapt.Filter")
        ylim([-.1 1.1])

        legend
        xlabel('Time (s)')        
        
    end
end

%% Effort Levels Plotting

Trials=[1 2];
ExpRep=5;
HalfTest=ExpRep*2;
TimeRange=[5 25];
eff=[];
for iTest=1:length(TestFolders)
    TestLabel=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpLabels(6);
    EffortFESTable=S.(TestLabel).(ExpLabel).EffortFESTable;
    CompTrialNum=sum(EffortFESTable(:,end));
    TrialNum=length(EffortFESTable(:,end));
    FrameRange=S.(TestLabel).stim_freq*TimeRange;
    FrameInd=FrameRange(1):FrameRange(2);    
    
    StimOnlyInd=logical(EffortFESTable(1:HalfTest,3));
    VoliOnlyInd=logical(EffortFESTable(1:HalfTest,4));
    EDCCFESInd=logical(EffortFESTable(HalfTest+1:HalfTest*2,2));
    CCFESInd=not(EffortFESTable(HalfTest+1:HalfTest*2,2));
    EffortTrials(:,1)=EffortFESTable(StimOnlyInd,1);
    EffortTrials(:,2)=EffortFESTable(VoliOnlyInd,1);
    EffortTrials(:,3)=EffortFESTable(EDCCFESInd,1)+HalfTest;
    EffortTrials(:,4)=EffortFESTable(CCFESInd,1)+HalfTest;
    TypeLabels=["VoliOnly" "StimOnly" "EDCCFES" "CCFES"];
    for iType=1:4
        TrialsNums=EffortTrials(:,iType);
        TypeLabel=TypeLabels(iType);

        for iTrial=1:length(TrialsNums)
            TrialLabel=sprintf("Trial_%d",TrialsNums(iTrial));

            Time=S.(TestLabel).(ExpLabel).(TrialLabel).time();
            TimeInd=Time>=TimeRange(1) & Time <=TimeRange(2);
            Hands=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,4:5);
            PW=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,6);
            Target=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,21);
            StimCont=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,22);
            Effort=S.(TestLabel).(ExpLabel).(TrialLabel).data(TimeInd,23);
            
            eff=[eff; sum(Hands(1)) sum(Hands(2)) sum(PW) sum(StimCont) sum(Hands(1))-sum(StimCont) sum(Effort) TrialsNums(iTrial) TypeLabel TestLabel]            ;
% figure
            % plot(Time(TimeInd),Effort)
             % rows are frames, columns are samples in a frame
            % ign=9;  
            % FiltvEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltvEMG_f(ign:end,:);
            % FiltmWaves_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltmWaves_f(ign:end,:);
            % Trig_f=S.(TestLabel).(ExpLabel).(TrialLabel).Trig_f(ign:end,:);
            % BlankedEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).BlankedEMG_f(ign:end,:);
            % Dropped_f=S.(TestLabel).(ExpLabel).(TrialLabel).DroppedEMG_f(ign:end,:);
            % DroppedTrig_f=S.(TestLabel).(ExpLabel).(TrialLabel).DrpdTrig_f(ign:end,:);
            % PW_f=S.(TestLabel).(ExpLabel).(TrialLabel).PW_f(ign:end,:);
            % 
            % Time_f=S.(TestLabel).(ExpLabel).(TrialLabel).Time_f(ign:end,:);                
            % Tm_f=S.(TestLabel).(ExpLabel).(TrialLabel).Tm_f(ign:end,:);
            % Time_f=Tm_f;
            % [r,c]=size(Time_f);
            % 
            % FiltvEMG_f=reshape(FiltvEMG_f',1,[]);
            % Time_f=reshape(Time_f',1,[]);
            % FiltmWaves_f=reshape(FiltmWaves_f',1,[]);
            % Trig_f=reshape(Trig_f',1,[]);
            % BlankedEMG_f=reshape(BlankedEMG_f',1,[]);                
            % Dropped_f=reshape(Dropped_f',1,[]);
            % DroppedTrig_f=reshape(DroppedTrig_f',1,[]);
            % PW_f=reshape(PW_f',1,[]);
            % 

        end     
    end
end


% efftable=array2table(eff,'VariableNames',["Paretic" "NonParetic" "PW" "StimLaw" "Error_tr" "Effort" "TrialNum" "TrialType" "Test"]);

% writetable(efftable,"effortFES2.csv")


%% smooth controller

% Define parameters
theta_np = 1;  % Example coefficient
E = 5;         % Example energy value
c = 2;         % Threshold for e_tr
h = 3;         % Threshold for Error_hand
k1 = 5;        % Smoothness factor for e_tr
k2 = 5;        % Smoothness factor for Error_hand

% Generate sample values for e_tr and Error_hand
e_tr_values = linspace(-5, 5, 100);
Error_hand_values = linspace(0, 5, 100);

% Store results
S_values = zeros(length(e_tr_values), length(Error_hand_values));

% Compute S_i for different values
for i = 1:length(e_tr_values)
    for j = 1:length(Error_hand_values)
        S_values(i, j) = smooth_S(theta_np, E, e_tr_values(i), c, Error_hand_values(j), h, k1, k2);
    end
end

% Plot the results
figure;
imagesc(Error_hand_values, e_tr_values, S_values);
colorbar;
xlabel('Error_{hand}');
ylabel('e_{tr}');
title('Smooth Transition of S_i');

% Define the function smooth_S
function S_i = smooth_S(theta_np, E, e_tr, c, Error_hand, h, k1, k2)
    % Smooth transition function for e_tr
    f_e = (1 / (1 + exp(-k1 * (e_tr - c)))) * min(1, e_tr / c);
    % f_e = (1 / (1 + exp(-k1 * e_tr))) * (1 - exp(-c * e_tr));
    
    % Ensure f_e is non-negative
    % f_e = max(f_e, 0);

    % Smooth transition function for Error_hand
    g_E = 1 / (1 + exp(-k2 * (Error_hand - h)));
    
    % Compute smooth S_i
    S_i = theta_np * E * f_e * g_E;
end


