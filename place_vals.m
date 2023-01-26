function array=place_vals(vals,ind,array)
% This function places values to an array while moving it 
% 1-values to be placed
% 2-indice of the values in array
% 3-array 


for i=1:length(ind)
    array(ind(i)+1:end+1)=array(ind(i):end);
    array(ind(i))=vals(i);
end


end