function S=tidy_data(FolderName,save)
%% This script is for tidying the data after a test with a participant
% Enter a string or char of the folder name where raw data is located. raw file
% name must be "expsave.mat"
% This will reorganize the data for further analyses 
% Only one file at a time whose name is expsave.mat


FolderName=string(FolderName);

FileName='expsave';    
ExpStruct=sprintf("%s_test",FolderName);
AnaStruct=sprintf('%s_ana',FolderName);

str=sprintf('%s/%s.mat',FolderName,FileName);
temp=load (str);
% S.(ExpStruct).ExpPar.Amp=temp.handles.Amp;
% S.(ExpStruct).ExpPar.MaxPW=temp.handles.MaxPW;
S.(ExpStruct).ExpPar.FreqList=temp.handles.freq_list;
% S.(ExpStruct).ExpPar.CalibMatrix=temp.handles.CalibMatrix;
S.(ExpStruct).ExpPar.FolderName=temp.handles.FolderName;
S.(ExpStruct).ExpPar.sample_t=temp.handles.sample_t;
S.(ExpStruct).ExpPar.fs=1/temp.handles.sample_t;
S.(ExpStruct).ExpPar.TrialsFreq=temp.handles.ExpTrials.TrialsFreq;
S.(ExpStruct).ExpPar.stim_freq=temp.handles.freq_list(temp.handles.ExpTrials.TrialsFreq(1));

TestExpNames=["MVCTrials","RCCurveTrials","CustomTrials","ExpTrials","FatigueTrials", "RCRampTrials"];
AnaExpNames=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials", "RCRampTrials"];

ExpID=field_exist(temp.handles,TestExpNames);

for iField=1:length(ExpID)
    AnaExpName=AnaExpNames(iField);
    TestExpName=TestExpNames(iField);
    if ExpID(iField)
        S.(ExpStruct).(AnaExpName)=temp.handles.(TestExpName); % if Exp exist then bring it to S structure
       
        if ~isfield(S.(ExpStruct).(AnaExpName), 'NumofTrials') % if NumofTrials does not exist then define it 
            S.(ExpStruct).(AnaExpName).NumofTrials=S.(ExpStruct).(AnaExpName).iTrial-1;
        end
        
        if S.(ExpStruct).(AnaExpName).NumofTrials==0 % if NumofTrials equal to 0 then Exp is never run thus ExpRun is false otherwise true
            S.(ExpStruct).(AnaExpName).ExpRun=false;
        else
            S.(ExpStruct).(AnaExpName).ExpRun=true;
        end
    else
        S.(ExpStruct).(AnaExpName).ExpRun=false;
        S.(ExpStruct).(AnaExpName).NumofTrials=0;
        S.(ExpStruct).(AnaExpName).iTrial=1;
    end
    
    ExpRuns(iField)=S.(ExpStruct).(AnaExpName).ExpRun;
end
S.(ExpStruct).ExpRuns=ExpRuns;


if isfield(temp.handles, 'SbjInfo')
    S.(ExpStruct).ExpPar.SbjInfo=temp.handles.SbjInfo;
end

if isfield(temp.handles, 'AmpGain')
    S.(ExpStruct).ExpPar.AmpGain=str2double(temp.handles.AmpGain);
end

if isfield(temp.handles, 'model')
    S.(ExpStruct).ExpPar.model=temp.handles.model;
else 
    S.(ExpStruct).ExpPar.model='EMGModel_test2';
end

if isfield(temp.handles, 'EffortType')
    
    S.(ExpStruct).ExpPar.EffortType=temp.handles.EffortType;
else
    S.(ExpStruct).ExpPar.EffortType="Force"; 
end

if isfield(temp.handles, 'EffortTypeButton')
    S.(ExpStruct).ExpPar.EffortTypeButton=temp.handles.EffortTypeButton;
end

if isfield(temp.handles, 'ExpLabels')
    S.(ExpStruct).ExpPar.ExpLabels=temp.handles.ExpLabels;

end

%
% Voli and Stim MVC 
%

AnaExpName= AnaExpNames(4); %% Occ Trials
% Only the test at feb29 has a non def stim vec

if FolderName=="feb29_24"
    DefStimMVCVec=[0 10 12 15 ];
else
    DefStimMVCVec=[0 10 20 30 ];
end
DefVoliMVCVec=[10 20 30 40];

if isfield(temp.handles.ExpTrials, 'StimMVCVec')
    if any(FolderName==["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25"])
        S.(ExpStruct).(AnaExpName).StimMVCVec=DefStimMVCVec;
    else
        S.(ExpStruct).(AnaExpName).StimMVCVec=temp.handles.ExpTrials.StimMVCVec;
    end
else 
    S.(ExpStruct).(AnaExpName).StimMVCVec=DefStimMVCVec;
end

if isfield(temp.handles.ExpTrials, 'VoliMVCVec')
    S.(ExpStruct).(AnaExpName).VoliMVCVec=temp.handles.ExpTrials.VoliMVCVec;
else 
    S.(ExpStruct).(AnaExpName).VoliMVCVec=DefVoliMVCVec;
end


%
% Labeling
%

if isfield(temp.handles, 'DataInd')
    S.(ExpStruct).ExpPar.DataInd=temp.handles.DataInd;

else 
    S.(ExpStruct).ExpPar.DataInd=["EMG","Trigger","Force","PW","Time", "Torque"];

end


S.(ExpStruct).ExpPar.DataIndTable=array2table(1:length(S.(ExpStruct).ExpPar.DataInd),'VariableNames',S.(ExpStruct).ExpPar.DataInd);
S.(ExpStruct).ExpPar.ExpLabels=AnaExpNames;
S.(ExpStruct).ExpPar.ExpLabelsTable=table(1,2,3,4,5,6,'VariableNames',AnaExpNames);
    
% 
% if isfield(temp.handles, 'ExpLabels')
%     S.(ExpStruct).ExpPar.ExpLabels=temp.handles.ExpLabels;
%     S.(ExpStruct).ExpPar.ExpLabelsTable=table(1,2,3,4,5,6,'VariableNames',S.(ExpStruct).ExpPar.ExpLabels);
% 
% else 
%     S.(ExpStruct).ExpPar.ExpLabels= ["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials"];
%     S.(ExpStruct).ExpPar.ExpLabelsTable=table(1,2,3,4,5,'VariableNames',S.(ExpStruct).ExpPar.ExpLabels);
% end
% 
% if isfield(temp.handles, 'RCRampTrials')
%     S.(ExpStruct).ExpPar.ExpLabels(end+1)="RCRampTrials";
%     S.(ExpStruct).ExpPar.ExpLabelsTable=table(1,2,3,4,5,6,'VariableNames',S.(ExpStruct).ExpPar.ExpLabels);
% 
% end
% ExpLabels=["MVCTrials","RCCurveTrials","CustomTrials","OccTrials","FatigueTrials" ];
% S.(ExpStruct).ExpPar.ExpLabels=ExpLabels;
% S.(ExpStruct).ExpPar.DataInd=table(1,2,3,4,5,6,'VariableNames',DataInd);

% 
% Num of Trials might not be needed anymore
%

ExptoAna=S.(ExpStruct).ExpPar.ExpLabels;

for iExp=1:length(ExptoAna)
    ExpLabel=ExptoAna(iExp);
    
    if ExpLabel=="OccTrials"
        NumofTrials=length(S.(ExpStruct).(ExpLabel).RepTableMat(:,6));
        S.(ExpStruct).(AnaExpName).ListedNumofTrials=NumofTrials;

    else
        NumofTrials=S.(ExpStruct).(ExpLabel).NumofTrials;
    end
    c=1;
    TrialNums=[];
    S.(ExpStruct).(ExpLabel).RedoTrials=TrialNums;
    for iTrial=1:NumofTrials
        RedoLabel=sprintf('RedoTrial_%d',iTrial);
        
        if isfield(S.(ExpStruct).(ExpLabel), RedoLabel)
            TrialNums(c)=iTrial;
            c=c+1;
        end
        S.(ExpStruct).(ExpLabel).RedoTrials=TrialNums;
        S.(ExpStruct).(ExpLabel).NumofRedos=length(TrialNums);
    end
    S.(ExpStruct).ExpPar.NumofTrialsTable{iExp,1}=ExpLabel;
    S.(ExpStruct).ExpPar.NumofTrialsTable{iExp,2}=S.(ExpStruct).(ExpLabel).NumofTrials;
end



%
% Here, check if certain variables exist and add if does not
%  experiments before feb27 need this part 
%


% StimRange for Exps


AmpGain=990;  % if AmpGain for nondefined
% S.(ExpStruct).ExpPar.AmpGain=AmpGain;
num_of_dropped=1;

% 1- AmpGain

if ~isfield(S.(ExpStruct).ExpPar, 'AmpGain')
    S.(ExpStruct).ExpPar.AmpGain=AmpGain;
end

% 2- dropped, num_of_dropped, of frames
dropped=0;
if ~isfield(S.(ExpStruct).RCCurveTrials, 'num_of_dropped')
    
    S.(ExpStruct).RCCurveTrials.num_of_dropped=num_of_dropped;
    S.(ExpStruct).RCCurveTrials.dropped=dropped;

end

dropped=0;
if ~isfield(S.(ExpStruct).CustomTrials, 'num_of_dropped')
    
    S.(ExpStruct).CustomTrials.num_of_dropped=num_of_dropped;
    S.(ExpStruct).CustomTrials.dropped=dropped;

end

dropped=1;
if ~isfield(S.(ExpStruct).OccTrials, 'num_of_dropped')
    
    S.(ExpStruct).OccTrials.num_of_dropped=num_of_dropped;
    S.(ExpStruct).OccTrials.dropped=dropped;

end

dropped=0;
if ~isfield(S.(ExpStruct).FatigueTrials, 'num_of_dropped')

    S.(ExpStruct).FatigueTrials.num_of_dropped=num_of_dropped;
    S.(ExpStruct).FatigueTrials.dropped=dropped;

end
    
% 3- stim_freq, FreqList

stim_freq=35;
if ~isfield(S.(ExpStruct).ExpPar, 'stim_freq')
    
    S.(ExpStruct).ExpPar.stim_freq=stim_freq;

end

FreqList=[ 35; 20; 50];
if ~isfield(S.(ExpStruct).ExpPar, 'FreqList')
    
    S.(ExpStruct).ExpPar.FreqList=FreqList;

end


if ~isfield(S.(ExpStruct).ExpPar, 'FreqList')
    
    S.(ExpStruct).ExpPar.FreqList=FreqList;

end

% change TrialsPW field name
if isfield(S.(ExpStruct).RCCurveTrials, 'PWTrials')
    
    TrialsPW=S.(ExpStruct).RCCurveTrials.PWTrials;
    S.(ExpStruct).RCCurveTrials.TrialsPW=TrialsPW;

end
   


% Stim Range for trials : might change for future trials


S.(ExpStruct).RCCurveTrials.StimRange=[5 10];
S.(ExpStruct).FatigueTrials.StimRange=[5 35];


%
% Consistent data-length
%

% RCTurnOffTimes= table([ 10 ],[10],[ 10],[8],[10],...
%     'VariableNames',["dec5","nov28_2","nov27","nov8","nov7"]);


mrg=[1.8 2 1];

sample_t=S.(ExpStruct).ExpPar.sample_t;
fs=1/sample_t;
TurnOffInd=4;


if ~isfield(S.(ExpStruct).OccTrials, 'StimProfile')
    S.(ExpStruct).OccTrials.StimProfile=15; % Stim turn off time 
end

if ~isfield(S.(ExpStruct).OccTrials, 'StimRange')
    S.(ExpStruct).OccTrials.StimRange=[S.(ExpStruct).OccTrials.StimProfile-10 S.(ExpStruct).OccTrials.StimProfile];
end

if ~isfield(S.(ExpStruct).OccTrials, 'StimConstantRange')
    S.(ExpStruct).OccTrials.StimConstantRange=[S.(ExpStruct).OccTrials.StimProfile-5 S.(ExpStruct).OccTrials.StimProfile];
end


S.(ExpStruct).FatigueTrials.StimProfile=35;
S.(ExpStruct).RCCurveTrials.StimProfile=10; % Stim turn off time 

Lbl=S.(ExpStruct).ExpPar.ExpLabels([2,4,5]);

for iExp=1:length(Lbl)
    ExpLabel=Lbl{iExp};

    T_turnoff=S.(ExpStruct).(ExpLabel).StimProfile;
    lgt=round(T_turnoff*fs+fs*mrg(iExp));
    
    NumofTrials=S.(ExpStruct).(ExpLabel).NumofTrials;
    RedoTrials=S.(ExpStruct).(ExpLabel).RedoTrials;

    for iTrial=1:NumofTrials
        TrialLabel=sprintf('Trial_%d',iTrial);

        % if ~isfield(S.(ExpStruct).OccTrials.(TrialLabel), 'StimProfile')
        %     S.(ExpStruct).OccTrials.(TrialLabel).StimProfile=15; % Stim turn off time 
        % end
        % 
        % if ~isfield(S.(ExpStruct).OccTrials.(TrialLabel), 'StimRange')
        %     S.(ExpStruct).OccTrials.(TrialLabel).StimRange=[S.(ExpStruct).OccTrials.(TrialLabel).StimProfile-10 S.(ExpStruct).OccTrials.(TrialLabel).StimProfile];
        % end
        % 
        % if ~isfield(S.(ExpStruct).OccTrials, 'StimConstantRange')
        %     S.(ExpStruct).OccTrials.(TrialLabel).StimConstantRange=[S.(ExpStruct).OccTrials.(TrialLabel).StimProfile-5 S.(ExpStruct).OccTrials.(TrialLabel).StimProfile];
        % end
        
        [lg_dt,~]=size(S.(ExpStruct).(ExpLabel).(TrialLabel).data);
        if  lg_dt>=lgt
            S.(ExpStruct).(ExpLabel).(TrialLabel).data(lgt:end,:)=[];
        else
            S.(ExpStruct).(ExpLabel).(TrialLabel).data(end+1:lgt,:)=0;
        end
    end 
    
    for iTrial=1:length(RedoTrials)
        RedoLabel=sprintf('RedoTrial_%d',RedoTrials(iTrial));
        
        [lg_dt,~]=size(S.(ExpStruct).(ExpLabel).(RedoLabel).data);
        if  lg_dt>=lgt
            S.(ExpStruct).(ExpLabel).(RedoLabel).data(lgt:end,:)=[];
        else
            S.(ExpStruct).(ExpLabel).(RedoLabel).data(end+1:lgt,:)=0;
        end    
    end
    
    S.(ExpStruct).ExpPar.lgt(iExp)=lgt;
end

S.(ExpStruct).ExpPar.fs_calc=fs;
S.(ExpStruct).ExpPar.sample_t_calc=sample_t;

ExpLabels=S.(ExpStruct).ExpPar.ExpLabels;  

ExpShortCuts=["MVC","RC","Cus","Occ","Fat","Ramp","Calib"];
ExpShortCutsTableVal=[ExpLabels; 1:length(ExpLabels)];
ExpShortCutsVarible=ExpShortCuts(1:length(ExpLabels));

S.(ExpStruct).ExpPar.ExpTable=array2table(ExpShortCutsTableVal,'VariableNames',ExpShortCutsVarible);

%
% Saving as struc
%
% 
% if nargin==2 & save==true
%     str=sprintf('%s/%s',FolderName,ExpStruct)
%     save(str,'-struct','S',ExpStruct)
%     str=sprintf("%s.mat is created",ExpStruct);
%     disp(str)
% end

end
