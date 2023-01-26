function []=save_test(TestFolders,S)
% This function saves the analysis results to their corresponding test folders
% enter folder names without the extension (such as _test )
for iTest=1:length(TestFolders)
    TestStruct=sprintf("%s_test",TestFolders{iTest});
    AnaStruct=sprintf("%s_ana",TestFolders{iTest});
    
    str=sprintf('%s/%s',TestFolders{iTest},AnaStruct);
    save(str,'-struct','S',TestStruct,AnaStruct)
end
end