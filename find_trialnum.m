function IndTrials=find_trialnum(vMVC,sMVC,RepTableMat)
% This function will find the occ trial numbers for a pair of given voli and stim MVC percent values
% vMVC: voli MVClevels such as 10, 20 etc
% sMVC: stim MVClevels such as 10, 20 etc
% RepTableMat: Occlusion experiment table found in Occ trials field
% IndTrials: array with trial numbers 
%----This needs revising 


    [VoliMVCTrials]=RepTableMat(:,4); 
    [StimMVCTrials]=RepTableMat(:,5); 
    % vMVC=VoliMVCLevels(PlotVoli);
    % sMVC=StimMVCLevels(PlotStim);
    IndV=find(VoliMVCTrials==vMVC);
    IndS=find(StimMVCTrials==sMVC); 
    
    if isempty (IndV)
        error('Volitional MVC value is not found ')
        
    end
    if isempty (IndV)
        error('Stim. MVC value is not found ')
        
    end

    for iV=1:length(IndV)
        x(IndS==IndV(iV))=1;
        IndTrials=IndS(x==1);
    end


end