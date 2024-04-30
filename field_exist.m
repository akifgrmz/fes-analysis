function res=field_exist(str,fieldnamess)

flnames=fieldnames(str); % fields to be checked 
fieldnamess=string(fieldnamess); % fieldnames that being looked for
res=logical(length(fieldnamess));

for iField=1:length(fieldnamess)
    
    temp=(flnames==fieldnamess(iField));
    res(iField)=any(temp);

end

