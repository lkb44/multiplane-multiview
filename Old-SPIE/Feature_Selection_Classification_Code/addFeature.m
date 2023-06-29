% Code for adding a single shape feature at a time
% Code written by Charlems

function feature_matrix_case = addFeature(imported_data, feature_matrix_case)

if ~contains(imported_data,'[')
   feature_matrix_case = [feature_matrix_case str2num(imported_data)];           
else
   tmp = strrep(imported_data,'[','');
   tmp = strrep(tmp,']','');
   tmp = strsplit(string(tmp),', ');

   for k=1:size(tmp,2) %move position by position into feature set [a, b, ... n]
       feature_matrix_case = [feature_matrix_case str2num(tmp{k})]; 
   end
end