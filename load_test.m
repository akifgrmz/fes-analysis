function M=load_test(FolderNames,FileName,StructNames)
% This function is designed to load the experiment data.
% Files must be in format defined previously.
% 1-enter an array of strings whose entries are folder names
% 2-enter an array of strings whose entries are file names
% 3-enter the name of the field when loading only spesific fields

%----> Modify so that it does load the tests already loaded


    if nargin ==2
        for iFile=1:length(FolderNames)
            ExpStruct=sprintf('%s',char(FileName{iFile}));
            str=sprintf('%s/%s',char(FolderNames{iFile}),ExpStruct);
            S=load (str);
            FldNames = fieldnames(S);
            for iField=1:length(FldNames)
                M.(FldNames{iField})=S.(FldNames{iField});
            end
        end

    elseif nargin ==3
        
        for iFile=1:length(FolderNames)
            ExpStruct=sprintf('%s',char(FileName{iFile}));
            str=sprintf('%s/%s',char(FolderNames{iFile}),ExpStruct);
            S=load (str);
            FldNames = fieldnames(S);
            for iStruct=1:length(StructNames)
                StructLabel=StructNames{iStruct};
                M.(ExpStruct).(StructLabel)=S.(ExpStruct).(StructLabel);
            end
        end

    elseif nargin == 1

        error('Enter a cell of strings that has at least length of 1 ')
    
    elseif nargin == 0
        
        for iFile=1:length(FolderNames)
            FileName=sprintf('%s_ana',char(FolderNames{iFile}));
            str=sprintf('%s/%s',char(FolderNames{iFile}),FileName);
            S=load (str);
            FldNames = fieldnames(S);
            for iField=1:length(FldNames)
                M.(FldNames{iField})=S.(FldNames{iField});
            end
        end
    end
end