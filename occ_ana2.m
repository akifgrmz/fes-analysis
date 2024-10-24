%% Occlusion Analysis
% This file contains all the occlusion related analyses
%% Data Inject 
clc
clear all
% TestFolders=["jan7" "jan11" "jan12" "apr20" "may19" "oct11" "oct18" "oct25"];
% TestFolders=["jan7" "jan11" "jan12" "apr20" ];
TestFolders=[  "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=[  "feb28_24" "mar20_24"];
TestFolders=[ "feb28_24" "feb29_24" "mar20_24"  ];
TestFolders=["jan7" "jan11" "jan12" "feb28_24"  ];
TestFolders=["jun20_24" "jul9_24" "jul21_24"];
TestFolders=["jan7" "jan11" "jan12"  ];
TestFolders=[ "jan7" "jan11" "jan12" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" ];
TestFolders=["feb28_24" "feb29_24"  "mar20_24" "aug19_24" ];
% TestFolders=["aug19_24" ];
TestFolders=["aug19_24" "aug22_24" "aug26_24"];

TestFolders=[ "jan7" "jan11" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" "jan12" "aug22_24" "aug26_24"];
TestFolders=[ "aug22_24" "aug26_24" "aug29_24"];
TestFolders=["jan7" "jan11" "jan12" "aug22_24" "sep3_24" "sep4_24" "sep6_24" "oct17_24" "oct18_24"];

TestFolders=["sep6_24" "oct17_24" "oct18_24"];

NumofTests=length(TestFolders);
DroppedFrameFilts=strings(1,NumofTests)+"GS";
NoStimFilts=strings(1,NumofTests)+"GS";
MAV_MAXMethods=strings(1,NumofTests)+"Real";

% AvgOcclusionTests=["jan7" "jan11" "jan12";
%      "aug22_24" "sep3_24" "sep4_24" ];
AvgOcclusionTests=["sep6_24" "oct17_24" "oct18_24"];

textt="%s_ana-oct19-";
TestFiles=compose(textt,TestFolders');
S = load_test(TestFolders,TestFiles);

Tau_table= readtable('tau_estimates3.csv');
Tau_stats_table= readtable('tau_est_stats.csv');

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

%% Effort Estimation
% Average Occlusion estimation
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

writetable( Effort_Est,'occ_est_error4.csv')

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

    sMVCzero=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('sMVC'); 
    MAVMean=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('MAV_Mean_Reps');
    vMVCLevs=S.(AnaLabel).(ExpLabel).MAV_Mean_table.('vMVC');
    
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

%% Plotting Occ Trials with Effort, Occlusion Estimation

close all
clc
PlotVoli=4;
PlotStim=4;
FiltLabel="GS";
iType=1;
TypeLabel=OccType(iType);
TimeRangeSNR=[6 11];
DroppedRange=2:6;
for iTest=1:length(TestFolders)
    AnaLabel=sprintf('%s_ana',TestFolders(iTest));
    TestLabel=sprintf('%s_test',TestFolders(iTest));
    EffortType=S.(TestLabel).ExpPar.EffortType;
    EffortLabel=sprintf("Filt_%s",EffortType);
    
    stim_freq=S.(TestLabel).ExpPar.stim_freq;
    FrameInd=stim_freq*TimeRange(1): stim_freq*TimeRange(2);
    FrameIndSNR=stim_freq*TimeRangeSNR(1): stim_freq*TimeRangeSNR(2);
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
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
    ttl=sprintf('Test: %s, sMVC: %d, vMVC: %d',TestFolders(iTest),sMVC,vMVC);
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('MVC');
    MAV_MAX=S.(AnaLabel).(ExpLabel).MAV_MAX_theo;
    MVC=S.(TestLabel).(ExpLabel).MVC;

    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).('Occ');
    stim_freq=S.(TestLabel).ExpPar.stim_freq;

    for iTrial=1:length(IndTrials)-2
        
        TrialLabel=sprintf('Trial_%d',IndTrials(iTrial));
        FiltForce=S.(AnaLabel).(ExpLabel).(TrialLabel).data.(EffortLabel);
        Time=S.(AnaLabel).(ExpLabel).(TrialLabel).data.('Time');        
        VoliLevel=RepTableMat(IndTrials(iTrial),3);
        Target=S.(AnaLabel).(ExpLabel).Target*RepTableMat(IndTrials(iTrial),1);
        DroppedMAV=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFeat(DroppedRange,:).('MAV');
        DroppedInd=S.(AnaLabel).(ExpLabel).(TrialLabel).DroppedFrames(DroppedRange)/stim_freq;

        E_f=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameIndSNR,:).('Effort_f');
        E_e_Amp(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_e_Amp');
        E_fv_prime(:,iTrial)=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).(TypeLabel).EffortEst(FrameInd,:).('Effort_Fv_prime');

        SNR= 20*log((RepTableMat(IndTrials(iTrial),4)+RepTableMat(IndTrials(iTrial),5))/std(E_f));

        figure(1)
        subplot(length(TestFolders),1,iTest)
        plot(Time,FiltForce,'DisplayName',sprintf('Force(N) (Trial: %d)',IndTrials(iTrial)))
        hold on
        plot(Time(1:end-1),Target,'DisplayName','Target Line')
        plot([10 15],[VoliLevel VoliLevel],'--','DisplayName','Voli Comp.')
        plot(FrameInd/stim_freq,E_e_Amp,'DisplayName','E_e')
        plot(FrameInd/stim_freq,E_fv_prime/100*MVC,'DisplayName','F_vprime')
        text(15.2,RepTableMat(IndTrials(iTrial),1)*12/10,sprintf('SNR:%.2f dB',SNR))
        ylim([-RepTableMat(IndTrials(iTrial),1)/10 RepTableMat(IndTrials(iTrial),1)*15/10])
    end
    
    title(ttl)
    legend('Location','NorthWest')
    grid on
end

%% Plotting for Conf Paper:EMBC
close all
clc
iRep=3;
TestFolders="jan7";
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
    FrameInd=stim_freq*TimeRange(1): stim_freq*TimeRange(2);

    ExpLabel=S.(TestLabel).ExpPar.ExpTable.(lbl);
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
TestFolders=["jun20_24" "jul9_24" ];
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
    
    ExpLabel=S.(TestLabel).ExpPar.ExpTable(1,:).(lbl);
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
    ExpLabel=S.(TestLabel).ExpPar.ExpTable.('Occ');
    
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
    







