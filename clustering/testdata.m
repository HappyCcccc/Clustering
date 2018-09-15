
%global evalCount;

%  Y1 = [randn(100,2)*0.75+ones(100,2);
%      randn(100,2)*0.5-ones(100,2)];
%   
%   for i=1:200
%       Z(i,1) = (Y1(i,1)-min(Y1(:,1)))/(max(Y1(:,1))-min(Y1(:,1)));
%       Z(i,2) = (Y1(i,2)-min(Y1(:,2)))/(max(Y1(:,2))-min(Y1(:,2)));
%   end

Z = rand(100,2);

%mesh(X,Y,cost);
opts = statset('Display','final');
[idx,D,sumd,D1] = kmeans(Z,2,'Distance','sqeuclidean',...
    'Replicates',1,'Options',opts);
%[D(:,1),D(:,2)]

% figure;
% plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',12)
% plot(D(:,1),D(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off


opts = statset('Display','final');
[idx,C] = kmeans(Z,2,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);
%[C(:,1),C(:,2)]

% figure;
% plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',12)
% hold on
% plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3) 
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off



