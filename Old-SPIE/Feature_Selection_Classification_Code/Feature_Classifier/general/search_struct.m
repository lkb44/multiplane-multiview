function indx_list=search_struct(structure,field,search)
indx_list=[];
for i=1:length(structure)
    if contains(getfield(structure(i),field),search,'IgnoreCase',true)
        indx_list(end+1)=i;
    end
end