
global Z
%0.6111
a = [0.7222 0.8333];
b = [0.1667 0.8333];
c = [0.1667 0.1667];

[a; b;c]
opts = statset('Display','final');
%[idx,C2] = kmeans(Z,3,'Distance','sqeuclidean','Replicates',1,'Options',opts);
[idx,C2] = kmeans(Z,3,'Replicates',1,'Options',opts);

[idx,C] = kmeans(Z,3,'Distance','sqeuclidean',...
    'Replicates',1,'Start',[a; b; c],'Options',opts);

C = [a b c]
figure;
plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',12)
hold on
plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',12)
hold on
plot(Z(idx==3,1),Z(idx==3,2),'g.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
 plot(C2(:,1),C2(:,2),'kx',...
    'color','blue', 'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
       'Location','SE')
title 'Cluster Assignments and Centroids'
hold on 
% plot(a(1),a(2),'r*','MarkerSize',15)
% hold on 
% plot(b(1),b(2),'r*','MarkerSize',15)
% hold on 
% plot(c(1),c(2),'r*','MarkerSize',15)
% 
hold off