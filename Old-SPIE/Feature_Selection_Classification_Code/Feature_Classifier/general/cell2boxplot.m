function cell2boxplot(data,f,t,xlab,ylab,leg)
alpha=length(data)
temp=[]
grp=[]
stack=[];
for i=1:alpha
[a,b]=size(data{i});
temp=i*ones(a,1);
grp=[grp;temp];
stackt=data{i}(:,f);
stack=[stack;stackt];
end


figure 
hold on
bp=boxplot(stack,grp)
% set(bp,'FontSize', '20')
set(bp,'LineWidth',2)
if nargin>2
    
bp=boxplot(stack,grp,'Labels',leg)
% set(bp,'FontSize', '20')
set(bp,'LineWidth',2)
title(t);
xlabel(xlab);
ylabel(ylab);
end
set(gca,'linewidth',2,'fontsize',24)
hold off
end