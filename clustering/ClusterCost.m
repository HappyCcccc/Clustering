
global Z



%global evalCount;

%  Y = [randn(10,2)*0.75+ones(10,2);
%      randn(10,2)*0.5-ones(10,2)];
% 
% for i=1:20
% Z(i,1) = (Y(i,1)-min(Y(:,1)))/(max(Y(:,1))-min(Y(:,1)));
% Z(i,2) = (Y(i,2)-min(Y(:,2)))/(max(Y(:,2))-min(Y(:,2)));
% end

k=10;
d=2;
w = 0.2*rand(1,k*d);
[X,Y] = meshgrid(0:0.01:1,0:0.01:1);
Z = randn(200,d);
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
