%% stats_dropped2

%% Stats with dropped frames
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders(iTest));
end

S = load_test(TestFolders,TestFiles);

%%
clc
MarginFromDropped=5;  % frames
TestStruct=sprintf("%s_test",TestFolders{iTest});

VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30];
stim_freq=S.(TestStruct).ExpPar.stim_freq;
TimeRange=[5 15];
FrameRange=TimeRange*stim_freq;  
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

% Normalization coeffs for the occ trials
sMVC=0;
vMVC=30;
MeanRange=[11 14];
MeanRangeInd=MeanRange*stim_freq;
for iTest=1:length(TestFolders)
    AnaStruct=sprintf("%s_ana",TestFolders(iTest));
    TestStruct=sprintf("%s_test",TestFolders(iTest));
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);
    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;

    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat) ;
    
    for iRep=1:length(IndTrials)
        RepLabel=sprintf("Trial_%d",IndTrials(iRep));
        
        MAV(:,iRep)=mean(S.(AnaStruct).(ExpLabel).(RepLabel)...
            .Unfilt.Feats.('MAV_vEMG')(MeanRangeInd(1):MeanRangeInd(2)));
        Voli_NormCoeff(:,iTest)=mean(mean(MAV));
        S.(AnaStruct).(ExpLabel).Voli_NormCoeff=Voli_NormCoeff(:,iTest);
    end 
end


FeatLabels=string(S.(AnaStruct).AnaPar.FeatLabels);
for iFeat=1:1
    FeatLabel=FeatLabels(iFeat);
    
    for iTest=1:length(TestFolders)
        TestStruct=sprintf("%s_test",TestFolders{iTest});
        AnaStruct=sprintf("%s_ana",TestFolders{iTest});

        RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
        Voli_NormCoeff(:,iTest)=S.(AnaStruct).(ExpLabel).Voli_NormCoeff;

        for iVoli=1:length(VoliMVCLevels)
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
                        
                        DroppedFeatLabel=sprintf('%s_vEMG',FeatLabel);
                        DroppedFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel);
                        DroppedFrameNum=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;                        
                        
                        
                        Feat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel);
                        Target=mean(S.(AnaStruct).(ExpLabel).TargetFrames)';
                        FrameNum=[1:length(Feat)];

                        x_feat=Feat;
                        x_target= Target(setdiff([1:length(Feat)+length(DroppedFeat)]',DroppedFrameNum));   
                        x_norm=x_feat/Voli_NormCoeff(:,iTest);
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
                        g_order(:)=NaN;
                        

                        g=[g; g_filt' g_voli' g_stim' g_test' g_feattype' g_rep' g_trial' g_order' ];
                        x=[x ;x_feat x_norm x_target x_framenum ];
                        
                    end
                    
                    DroppedFeatLabel=sprintf('%s_vEMG',FeatLabel);
                    num_of_dropped=S.(TestStruct).(ExpLabel).num_of_dropped;
                    DroppedFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel);
                    DroppedFrameNum=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                    DroppedOrder=repmat([1:num_of_dropped],1,length(DroppedFeat)/num_of_dropped);

                    x_droppedfeat=DroppedFeat;
                    x_droppedtarget=Target(DroppedFrameNum);
                    x_droppednorm=x_droppedfeat/Voli_NormCoeff(:,iTest);
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
% % DirLabelCSV=sprintf('%s/%s_dropped.csv',TestFolders{iTest},TestFolders{iTest});
writetable( Dropped_stats, 'dropped_stats2.csv')
% save_test(TestFolders,S)



















