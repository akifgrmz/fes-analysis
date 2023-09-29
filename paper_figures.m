%% Effort Paper Figures 
% # Data Inject
clc
clear all
TestFolders=["jan7" "jan11" "jan12" "apr20" "may19"];

for iTest=1:length(TestFolders)
    TestFiles(iTest)=sprintf("%s_ana",TestFolders{iTest});
end

S = load_test(TestFolders,TestFiles);
%%

