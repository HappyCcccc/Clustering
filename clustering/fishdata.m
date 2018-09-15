global Z;
global k;
global d1;
global datan;

load fisheriris
Z = meas(:,3:4);

% figure;
% plot(X(:,1),X(:,2),'k*','MarkerSize',5);
% title 'Fisher''s Iris Data';
% xlabel 'Petal Lengths (cm)'; 
% ylabel 'Petal Widths (cm)';

opts = statset('Display','final');
[idx,C] = kmeans(Z,2,'Distance','sqeuclidean',...
    'Replicates',1,'Options',opts);

opts = statset('Display','final');
[idx,C] = kmeans(Z,2,'Distance','sqeuclidean',...
    'Replicates',10,'Options',opts);