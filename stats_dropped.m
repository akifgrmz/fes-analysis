%% Stats with dropped frames
clc
clear all
TestFolders=["jan7" "jan11" "jan12" ];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders(iTest));
end

S = load_test(TestFolders,TestFiles);

%%
clc
MarginFromDropped=5;  % frames
TestStruct=sprintf("%s_test",TestFolders{iTest});
PlotRange1=[ 5.4 10]; 
PlotRange2=[ 10 15.2]; 
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0];

stim_freq=S.(TestStruct).ExpPar.FreqList(1);
PlotRangeFrames1=PlotRange1*stim_freq;  
PlotRangeFrames2=PlotRange2*stim_freq;
boolean DroppedFrameInd;
TrialColor={'r', 'b', 'k'};
exp_lbl='Occ';
VarNames=["Dropped" "Frame"  "Filt_Type" "Repeat"  "Trial"...
    "MVC_Voli" "MVC_Stim" "Feat" "Test" ];
NormVoliMVC=20;
g=strings(1,length(VarNames))+NaN;
x=ones(0,3);
x_frame=[];

% Normalization coeffs for the occ trials
sMVC=0;
vMVC=30;
MeanRange=[11 12];
MeanRangeInd=MeanRange*stim_freq;
for iTest=1:length(TestFolders)
    AnaStruct=sprintf("%s_ana",TestFolders(iTest));
    TestStruct=sprintf("%s_test",TestFolders(iTest));
    ExpLabel=S.(AnaStruct).AnaPar.ExpTable.(exp_lbl);
    RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;

    IndTrials=find_trialnum(vMVC,sMVC,RepTableMat) ;
    
    for iRep=1:length(IndTrials)
        RepLabel=sprintf("Trial_%d",IndTrials(iRep));
        
        MAV(:,iRep)=mean(S.(AnaStruct).(ExpLabel).(RepLabel).Unfilt.Feats.('MAV_vEMG')(MeanRangeInd(1):MeanRangeInd(2)));
        Voli_NormCoeff(:,iTest)=mean(mean(MAV));
        S.(AnaStruct).(ExpLabel).Voli_NormCoeff=Voli_NormCoeff(:,iTest);
    end 
end

for iFeat=1:1
    for iTest=1:length(TestFolders)
        TestStruct=sprintf("%s_test",TestFolders{iTest});
        AnaStruct=sprintf("%s_ana",TestFolders{iTest});

        FeatLabel=S.(AnaStruct).AnaPar.FeatLabels{iFeat};
        RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;
        Target=mean(S.(AnaStruct).(ExpLabel).TargetFrames)';
        Voli_NormCoeff(:,iTest)=S.(AnaStruct).(ExpLabel).Voli_NormCoeff;

        for iVoli=1:length(VoliMVCLevels)
            for iStim=1:length(StimMVCLevels)
                sMVC=StimMVCLevels(iStim);
                vMVC=VoliMVCLevels(iVoli);
                IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

                for iRep=1:length(IndTrials)
                    TrialLabel=sprintf('Trial_%d',IndTrials(iRep));
                    DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;
                    DroppedFeatLabel=sprintf('%s_vEMG',FeatLabel);

                    DroppedFrameInd1=PlotRangeFrames1(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames1(2);
                    DroppedFrameInd2=PlotRangeFrames2(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames2(2);
                    DroppedFrameInd=DroppedFrameInd1+DroppedFrameInd2;
                    DroppedFrameInd=logical(DroppedFrameInd);
                    DroppedFramesFeat=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel)(DroppedFrameInd);

                    if isempty(DroppedFrames)
                        FiltLabels=["Unfilt"];
                        DroppedFramesFeat=[];
                        DroppedFrames=
                    else
                    end

                    DroppedFlag=true(length(DroppedFrames),1);
                    for iFilt=1:length(FiltLabels)
                        FiltLabel=FiltLabels{iFilt};
                        vEMGLabel=sprintf('%s_vEMG',FeatLabel);

                        Feat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel);
                        FeatwithDropped=place_vals(DroppedFramesFeat,DroppedFrames,Feat);

                        for iDropped=1:length(DroppedFrames(DroppedFrameInd))
                            DroppedLabel=sprintf('Dropped_%d',iDropped);
                            DropLabel=sprintf('%d',iDropped);

                            DroppedInd=DroppedFrames(iDropped);
                            PlotInd=[-MarginFromDropped+DroppedInd:DroppedInd-1 DroppedInd+1:DroppedInd+MarginFromDropped];
%                             OnesInd=ones(length(PlotInd),1);

                            RelatedFrames=FeatwithDropped(PlotInd);
                            dp=[];
                            x_target=Target(PlotInd);
                            
                            if DroppedFlag(iDropped)
                                RelatedFrames(end+1)=DroppedFramesFeat(iDropped); 
                                dp=DroppedInd; 
                                x_target(end+1)=Target(DroppedInd);
                            end
                           
                            x_frame=RelatedFrames;
                            x_norm=x_frame/Voli_NormCoeff(:,iTest);

                            lg=length(x_frame);
                            g_filt=strings(1,lg);
                            g_voli=strings(1,lg);
                            g_stim=strings(1,lg);
                            g_test=strings(1,lg);
                            g_feat=strings(1,lg);
                            g_rep=strings(1,lg);
                            g_trial=strings(1,lg);
                            g_dropped=strings(1,lg);
                            g_frame=strings(1,lg);

                            g_voli(:)=sprintf("%d%%",(VoliMVCLevels(iVoli)));
                            g_stim(:)=sprintf("%d%%",(StimMVCLevels(iStim)));
                            g_rep(:)=sprintf("Rep_%d",iRep);
                            g_filt(:)=string(FiltLabel);
                            g_trial(:)=string(TrialLabel);
                            g_test(:)=string(TestFolders(iTest));
                            g_feat(:)=string(FeatLabel);
                            g_dropped(:)=string(DropLabel);
                            g_frame(:)=string([PlotInd dp ]);

                            if DroppedFlag(iDropped), g_filt(end)="Dropped"; DroppedFlag(iDropped)=false; end

                            g=[ g; g_dropped' g_frame' g_filt' g_rep' g_trial' g_voli' g_stim' g_feat' g_test' ];
                            x=[x ;x_frame x_norm x_target];
                        end
                    end
                end
            end
        end
    end
end


Dropped_stats=array2table([x g(2:end,:) ],'VariableNames',["Frame_Val" "Norm_Val" "Target" VarNames]);
S.(AnaStruct).(ExpLabel).Dropped_stats=Dropped_stats;
% DirLabelCSV=sprintf('%s/%s_dropped.csv',TestFolders{iTest},TestFolders{iTest});
writetable( Dropped_stats, 'dropped_stats.csv')

% save_test(TestFolders,S)














