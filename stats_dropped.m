%% Stats with dropped frames
clc
clear all
TestFolders=["jan7" "jan11" "jan12" ];

for iTest=1:length(TestFolders)
    TestFiles{iTest}=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);

%%
clc
MarginFromDropped=5;  % frames
TestStruct=sprintf("%s_test",TestFolders{iTest});
PlotRange1=[ 5.4 10]; 
PlotRange2=[ 10 15.2]; 
VoliMVCLevels=[10 20 30 40 ];
StimMVCLevels=[ 0 10 20 30];
stim_freq=S.(TestStruct).ExpPar.FreqList(1);
PlotRangeFrames1=PlotRange1*stim_freq;  
PlotRangeFrames2=PlotRange2*stim_freq;
boolean DroppedFrameInd;
TrialColor={'r', 'b', 'k'};
exp_lbl='Occ';
VarNames=["Dropped" "Frame"  "Filt_Type" "Repeat"  "Trial" "MVC_Voli" "MVC_Stim" "Feat" "Test" ];
g=strings(1,length(VarNames))+NaN;
x=[];
x_frame=[];

for iFeat=1:1
    for iTest=1:length(TestFolders)
        TestStruct=sprintf("%s_test",TestFolders{iTest});
        AnaStruct=sprintf("%s_ana",TestFolders{iTest});

        ExpLabel=string(S.(AnaStruct).AnaPar.ExpTable.(exp_lbl));
        FeatLabel=S.(AnaStruct).AnaPar.FeatLabels{iFeat};
        RepTableMat=S.(TestStruct).(ExpLabel).RepTableMat;

        for iVoli=1:length(VoliMVCLevels)
            for iStim=1:length(StimMVCLevels)
                sMVC=StimMVCLevels(iStim);
                vMVC=VoliMVCLevels(iVoli);
                IndTrials=find_trialnum(vMVC,sMVC,RepTableMat);

                for iRep=1:length(IndTrials)
                    TrialLabel=sprintf('Trial_%d',IndTrials(iRep));
                    DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                    TargetFrames=mean(S.(AnaStruct).(ExpLabel).TargetFrames);
                    FiltLabels=S.(AnaStruct).AnaPar.FiltLabels;

                    if isempty(DroppedFrames)
                        continue
                    end
                    DroppedFlag=true(length(DroppedFrames),1);
                    for iFilt=1:length(S.(AnaStruct).AnaPar.FiltLabels)
                        FiltLabel=S.(AnaStruct).AnaPar.FiltLabels{iFilt};
                        vEMGLabel=sprintf('%s_vEMG',FeatLabel);
                        DroppedFeatLabel=sprintf('%s_vEMG',FeatLabel);

                        DroppedFrames=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFrames;
                        DroppedFrameInd1=PlotRangeFrames1(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames1(2);
                        DroppedFrameInd2=PlotRangeFrames2(1)<= DroppedFrames & DroppedFrames<=PlotRangeFrames2(2);
                        DroppedFrameInd=DroppedFrameInd1+DroppedFrameInd2;
                        DroppedFrameInd=logical(DroppedFrameInd);

                        DroppedFramesFeat(iRep,:)=S.(AnaStruct).(ExpLabel).(TrialLabel).DroppedFeat.(DroppedFeatLabel)(DroppedFrameInd);

                        for iDropped=1:length(DroppedFrames(DroppedFrameInd))
                            DroppedLabel=sprintf('Dropped_%d',iDropped);
                            DropLabel=sprintf('Drop_%d',iDropped);

                            DroppedInd=DroppedFrames(iDropped);
                            PlotInd=[-MarginFromDropped+DroppedInd:DroppedInd-1 DroppedInd+1:DroppedInd+MarginFromDropped];
                            OnesInd=ones(length(PlotInd),1);

                            Feat=S.(AnaStruct).(ExpLabel).(TrialLabel).(FiltLabel).Feats.(vEMGLabel);
                            FeatwithDropped=place_vals(DroppedFramesFeat,DroppedFrames,Feat);
                            RelatedFrames=FeatwithDropped(PlotInd);
                            dp=[];
                            if DroppedFlag(iDropped), RelatedFrames(end+1)=DroppedFramesFeat(iRep,iDropped); dp=DroppedInd; end
                           
                            x_frame=RelatedFrames';
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
                            x=[x x_frame];
                        end
                    end
                end
            end
        end
    end
end
Dropped_stats=array2table([x' g(2:end,:) ],'VariableNames',["Frame_Val" VarNames]);
S.(AnaStruct).(ExpLabel).Dropped_stats=Dropped_stats;
% DirLabelCSV=sprintf('%s/%s_dropped.csv',TestFolders{iTest},TestFolders{iTest});
writetable( Dropped_stats, 'dropped_stats.csv')

%%


[p.(FiltLabel)(iVoli,iStim),t,stats]= anova1(y_frame,g_frame,'off');
results.(FiltLabel){iVoli,iStim}=multcompare(stats,'Display','off');
results.p.(FiltLabel)(iVoli,iStim)= results.(FiltLabel){iVoli,iStim}(3,6);















