%% Manual fixes and additions

% Run once for each experiment
manu_add;


%% Define Tests and Parameters 

% TestFolders=["mar20_24"];
% TestFolders=["jan7" "jan11" "jan12" "feb27" "mar7" "mar16" "apr20" "oct18" "oct25" "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=[ "apr20" "oct11" "oct18"];
% TestFolders=["oct25"];
% TestFolders=[ "feb28_24" "feb29_24" "mar18_24" "mar20_24"];
% TestFolders=["jan7" "jan11" "jan12"  ];
% TestFolders=["jun20_24" "jul9_24" "jul21_24"  ]
% TestFolders=[ "jan7" "jan11" "jan12" "jun20_24" "jul9_24" "jul21_24" "jul31_24" ];
TestFolders=["feb28_24" "feb29_24" "mar18_24"  "mar20_24" ];
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24"];
TestFolders=[ "jan7" "jan11" "jun20_24" "jul9_24" "jul21_24" "jul31_24" "aug19_24" "jan12" "aug22_24" "aug26_24"];
TestFolders=["aug29_24"]

%%  
TestFolders=[ "jan7" "jan11" "jan12" "aug22_24" "aug26_24"];

DroppedFrameFilts=strings(1,length(TestFolders))+"GS";
NoStimFilts=strings(1,length(TestFolders))+"Unfilt";
MAV_MAXMethods=strings(1,length(TestFolders))+"Fitted";

TauTestsForce=["jan7" "jan11" "jan12"]; % Average of 
TauTestsHand=["jan7" "jan11" "jan12"]; % Average of 


