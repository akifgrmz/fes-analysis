function []=save_test(Tests,S)
% Type: save_test(TestFolders,S) to save the analysis.
%This function saves the analysis results to their corresponding test
% folders 
% enter folder names without the extension (such as _test )


for iTest=1:length(Tests)
    TestStruct=sprintf("%s_test",Tests{iTest});
    AnaStruct=sprintf("%s_ana",Tests{iTest});

    str=sprintf('%s/%s',Tests{iTest},AnaStruct);
    save(str,'-struct','S',TestStruct,AnaStruct)
    fprintf("%s is Saved (%d/%d)\n",Tests{iTest},iTest,length(Tests))
end
    
end