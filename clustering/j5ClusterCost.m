
global Z



%global evalCount;

%  Y = [randn(10,2)*0.75+ones(10,2);
%      randn(10,2)*0.5-ones(10,2)];
% 
% for i=1:20
% Z(i,1) = (Y(i,1)-min(Y(:,1)))/(max(Y(:,1))-min(Y(:,1)));
% Z(i,2) = (Y(i,2)-min(Y(:,2)))/(max(Y(:,2))-min(Y(:,2)));
% end

k=5;
d=2;
w = 0.2*rand(1,k*d);
[X,Y] = meshgrid(0:0.01:1,0:0.01:1);
%Z = randn(20,d);
Z = zeros(100,d);
gap = 0.0;
iter = 0;
while (gap < 0.1)
for i=1:20
    iter = iter +1;
    tz = [-1 -1];
    while ((tz(1) < 0.0) || (tz(1) > 1.0) || (tz(2) < 0.0) || (tz(2) > 1.0))
      tz = mvnrnd([0.7 0.9], [0.01 0.005; 0.005 0.01]);
    end;
    Z(i,1) = tz(1);
    Z(i,2) = tz(2);
end;
for i=21:40
    tz = [-1 -1];
    while ((tz(1) < 0.0) || (tz(1) > 1.0) || (tz(2) < 0.0) || (tz(2) > 1.0))
      tz = mvnrnd([0.1 0.05], [0.01 0.005; 0.005 0.01]);
    end;
    Z(i,1) = tz(1);
    Z(i,2) = tz(2);
end;
for i=41:60
    tz = [-1 -1];
    while ((tz(1) < 0.0) || (tz(1) > 1.0) || (tz(2) < 0.0) || (tz(2) > 1.0))
      tz = mvnrnd([0.1 0.75], [0.01 0.005; 0.005 0.01]);
    end;
    Z(i,1) = tz(1);
    Z(i,2) = tz(2);
end;
for i=61:80
    tz = [-1 -1];
    while ((tz(1) < 0.0) || (tz(1) > 1.0) || (tz(2) < 0.0) || (tz(2) > 1.0))
      tz = mvnrnd([0.8 0.5], [0.01 -0.005; -0.005 0.01]);
    end;
    Z(i,1) = tz(1);
    Z(i,2) = tz(2);
end;
for i=81:100
    tz = [-1 -1];

    while ((tz(1) < 0.0) || (tz(1) > 1.0) || (tz(2) < 0.0) || (tz(2) > 1.0))
      tz = mvnrnd([0.5 0.5], [0.02 -0.01; -0.01 0.02]);
    end;
    Z(i,1) = tz(1);
    Z(i,2) = tz(2);
end;
opts = statset('Display','final');
[idx,C1, sumd1] = kmeans(Z,5,'Replicates',1,'Options',opts);
[idx,C2, sumd2] = kmeans(Z,5,'Replicates',100,'Options',opts);
Z
gap = sum(sumd1)-sum(sumd2)
'relative error'
gap/sum(sumd1)
end;
iter


figure;
plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',12)
hold on
plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',12)
hold on
plot(Z(idx==3,1),Z(idx==3,2),'g.','MarkerSize',12)
hold on
plot(Z(idx==4,1),Z(idx==4,2),'c.','MarkerSize',12)
hold on
plot(Z(idx==5,1),Z(idx==5,2),'m.','MarkerSize',12)
plot(C1(:,1),C1(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
       'Location','SE')
title 'Cluster Assignments and Centroids'
hold off


figure;
plot(Z(idx==1,1),Z(idx==1,2),'r.','MarkerSize',12)
hold on
plot(Z(idx==2,1),Z(idx==2,2),'b.','MarkerSize',12)
hold on
plot(Z(idx==3,1),Z(idx==3,2),'g.','MarkerSize',12)
hold on
plot(Z(idx==4,1),Z(idx==4,2),'c.','MarkerSize',12)
hold on
plot(Z(idx==5,1),Z(idx==5,2),'m.','MarkerSize',12)
plot(C2(:,1),C2(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Centroids',...
       'Location','SE')
title 'Cluster Assignments and Centroids'
hold off


%scatter(Z(:,1), Z(:,2))
return
pause(20)
for i=1:length(X)
    for j=1:length(Y)
        %(i,j) = (X(1,i)-0.5)^2 + (Y(j,1)-0.5)^2;
        alp = X(1,i);
        bet =  Y(j,1);
        w(k*d-1) = alp;
        w(k*d) = bet;
        cost(i,j) = cluster(w);
    end;
end;

mesh(X,Y,cost)
