function list= dir2(filepath)
dirlist=dir(filepath);
c=1;

for i=1:length(dirlist)
    if ~contains(dirlist(i).name,'.')&&~contains(dirlist(i).name,['.*'])&&~contains(dirlist(i).name,'.ini')
    list(c)=dirlist(i);
    c=c+1;
    end
end

end