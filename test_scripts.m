%% Various Test Scripts 

clear all
TestFolders="mar7";
% "jan11" "jan12" "feb27" "mar7" "mar16" "apr20"
AnaLabel=sprintf("%s_ana",TestFolders);

S = load_test(TestFolders,AnaLabel);
%%
FiltLabel="GS";
TrialNum=12;

scale=701;
offset=0.00001;
asfMultip=4.886;
TimeRange=[6 15];
TestLabel=sprintf('%s_test',TestFolders(1));
ExpLabel=S.(AnaLabel).AnaPar.ExpTable.('Occ');
stim_freq=S.(TestLabel).ExpPar.stim_freq;
FrameInd=ceil([stim_freq*TimeRange(1): stim_freq*TimeRange(2)]);
TrialLabel=sprintf("Trial_%d",TrialNum);
MAV_vEMG=S.(AnaLabel).(ExpLabel).(TrialLabel).(FiltLabel).Feats(FrameInd,:).('MAV_vEMG');
FrameInd2=ceil([stim_freq*10: stim_freq*15]);
MeanMAV=mean(MAV_vEMG(end-stim_freq*5:end));
asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11]*MeanMAV*asfMultip;

dir=0;
asf_index=1;

% for iFrame=2:length(MAV_vEMG)
%     
%     [Amp_MAV_vEMG(iFrame),dir,asf_index]=...
%         adap_filt_test(MAV_vEMG(iFrame-1:iFrame),scale,offset,asf,dir,asf_index);
% end

[Amp_MAV_vEMG,Clip_MAV_vEMG]=amp_modul(MAV_vEMG,scale,offset,asf); 

figure(1)
plot(Amp_MAV_vEMG)
hold on 
plot(MAV_vEMG*100)

% plot(MAV_vEMG)




function [out,dir,asf_index]=adap_filt_test(InSig,scale,offset,asf,dir,asf_index)
% Adaptive amplitude modulation for low pass filtering MAV signals

def_asf=[1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,5,5,6,7,8,9,9,10,11,11];


def_offset=1;
def_scale=0;

max_index=length(asf);
% ClippedSig=clip((InSig-offset)*scale,0,1);
% OutSig=ClippedSig*0;
OutSig(1)=InSig(1);
OutSig(2)=0;
%     delta=ClippedSig(2)-OutSig(1);
    delta=OutSig(2)-OutSig(1);
    new_dir=sign(delta);
    mag=abs(delta);

    if new_dir==0
        asf_index=1;
    elseif new_dir==dir
        asf_index=asf_index+1;

        if asf_index>max_index
            asf_index=max_index;
        end
    else
        asf_index=1;
    end
    
    dir=new_dir;
    step=min(mag,asf(asf_index))*dir;
    OutSig(2)=OutSig(1)+step;
    
%     OutSig=clip(OutSig,0,.15);
out=OutSig(2);

% function[clipped]=clip(input,low,high)
%     clipped=input*0;
%     for i=1:length(input)
%         clipped(i)=min(max(input(i),low),high);
%     end
% end

end