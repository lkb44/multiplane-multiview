%% PURPOSE OF THIS ANALYSIS: How are features/predictions holding up across
%slices for a particular patient?

%***ONLY MEANT TO BE CALLED WITHIN nFOLDCV_withFS_v3.m

p = tes{i}; %current patient
l = unique(testing_labels); %class label for that patient

if i<10
    t = sprintf('Iter=%i,Fold=0%i Testing Patient %i (label=%i)',j,i,pts(p),l);
else
    t = sprintf('Iter=%i,Fold=%i Testing Patient %i (label=%i)',j,i,pts(p),l);
end

figure('visible','off','units','normalized','position',[0 0 1 1])
% suptitle(t);

subplot(2,2,2);
plot(testing_set(:,1),'r--o','Markersize',12,'linewidth',2);xlabel('slices');ylabel('normalized feature value');xlim([0 length(testing_set(:,fea(1)))+1]);ylim([-3 3]);hold on;
plot(testing_set(:,2),'g--o','Markersize',12,'linewidth',2);xlabel('slices');ylabel('normalized feature value');xlim([0 length(testing_set(:,fea(1)))+1]);ylim([-3 3]);hold on;
plot(testing_set(:,3),'b--o','Markersize',12,'linewidth',2);xlabel('slices');ylabel('normalized feature value');xlim([0 length(testing_set(:,fea(1)))+1]);ylim([-3 3]);hold on;
title('Testing patient: Feature value along slices');set(gca,'Xtick',1:length(testing_set(:,fea(1))),'Ytick',-3:0.5:3);
legend(['TopFeat1--' int2str(fea(1))], ['TopFeat2--' int2str(fea(2))],['TopFeat3--' int2str(fea(3))]);

subplot(2,2,4);
plot(temp_stats.prediction,'k--o','Markersize',12,'linewidth',2); xlabel('slices');ylabel('classifier predictions');xlim([0 length(testing_set(:,fea(1)))+1]);
title('Testing Patient:Classifier Predictive Value along slices');set(gca,'Xtick',1:length(testing_set(:,fea(1))));ylim([0 1]);

% subplot(3,1,3);
% plot(temp_stats.prediction); xlabel('slices');ylabel('classifier predictions');
% title('Classifier Predictive Value along slices');


%% PURPOSE OF THIS ANALYSIS: What kind of classifier does it appear we will need?

% Visualize the partitions
[x,y,z] = meshgrid(linspace(-3,3,20),linspace(-3,3,20),linspace(-3,3,20));
x = x(:); y = y(:); z = z(:);
[C,err,P,logp,coeff] = classify([x y z],training_set(:,fea),training_labels,'linear');
K = coeff(1,2).const;
L = coeff(1,2).linear;
fx = @(x) K + x*L(1);
fy = @(y) K + y*L(2);
fz = @(z) K + z*L(3);

subplot(2,2,1);
scatter3(training_set((training_labels==1),fea(1)),training_set((training_labels==1),fea(2)),training_set((training_labels==1),fea(3)),'rs');hold on;
scatter3(training_set((training_labels==-1),fea(1)),training_set((training_labels==-1),fea(2)),training_set((training_labels==-1),fea(3)),'bd');
title('Training Data');
legend('Pos (pCR)','Neg (non-pCR');
xlabel('TopFeat1');ylabel('TopFeat2');zlabel('TopFeat3');
set(gca,'Xlim',[-3 3],'Ylim',[-3 3],'Zlim',[-3 3]);
view(-20,35);
hold off;


subplot(2,2,3);
scatter3(x(C==1),y(C==1),z(C==1),'r.');hold on;
scatter3(x(C==-1),y(C==-1),z(C==-1),'b.');
title('Training Model Space');
legend('Pos (pCR)','Neg (non-pCR');
xlabel('TopFeat1');ylabel('TopFeat2');zlabel('TopFeat3');
set(gca,'Xlim',[-3 3],'Ylim',[-3 3],'Zlim',[-3 3]);
% h = fplot3(fx,fy,fz);
% set(h,'Color','m','LineWidth',2)
view(-20,35);
hold off;

savefolder = 'G:\Team Drives\INVent\uh_rectal_radiology\temp_experiments\pcr_vs_nonpcr\results\temp_analysis\';

print([savefolder t],'-dpng');