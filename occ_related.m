%% Occlusion Analysis
% This file contains all the occlusion related analyses
%% Data Inject 
clear all
S=load_test;
%% Plotting the Occ Trials
close all
clc
lbl='Occ';
PlotVoli=1;
PlotStim=4;
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30 ];

for iTest=1:length(TestFolders)
    AnaStruct=sprintf('%s_ana',TestFolders(iTest));
    TestStruct=sprintf('%s_test',TestFolders(iTest));
    
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(lbl);
    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
    sMVC=StimMVCLevels(PlotStim);
    vMVC=VoliMVCLevels(PlotVoli);
    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

    StimOff=S.(TestStruct).(ExpLabel).TargetProfile2(4);     %this part might need change based new updates 
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


%% Time constants