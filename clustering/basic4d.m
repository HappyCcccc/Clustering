clear all; close all; clc;

global X;
global evalCount;

 Y = [randn(100,2)*0.75+ones(100,2);
     randn(100,2)*0.5-ones(100,2)];

for i=1:200
X(i,1) = (Y(i,1)-min(Y(:,1)))/(max(Y(:,1))-min(Y(:,1)));
X(i,2) = (Y(i,2)-min(Y(:,2)))/(max(Y(:,2))-min(Y(:,2)));
end


opts = statset('Display','final');
[idx,C] = kmeans(X,2,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);
[C(:,1),C(:,2)]

r = 1.5;  %smoothness parameter
evalCount = 0;
%z =0.0;
% Compute the value of the quadratic function.
tic
d = 4;
n = 10000;
Corner = zeros(1,d);% matrix for corner node
Width = zeros(1,d);% matrix for width and length
Center = zeros(1,d);% matrix for center node
area = ones(1,1); % matrix for retrangle area
f = zeros(1,1); % function value for center points
getrho = zeros(1,1);
a = zeros(1,d);
b = zeros(1,d);
V = zeros(1,n);
M = zeros(1,n);
x = zeros(1,n);
y = zeros(1,n);
z = zeros(1,n);
%A = zeros(1,d);


for i = 1:d
    Width(1,i)=1.0;
    Center(1,i)=0.5;
end

maxWidth = 0.0;
maxIndex = 1;
maxD = 0;

Mn = intmax('int64');
Vn = intmax('int64');
Center(maxIndex,:);

for i = 1:n
    
maxWidth = 0.0;
% split the retrangle
    %decide which axis to split on
for j = 1:d
        if(Width(maxIndex,j)>maxWidth)
            maxWidth = Width(maxIndex,j);
            maxD = j;
        end
end

% A =[];
% for j=1:d
%    if(Width(maxIndex,j)>0.9*maxWidth) 
%     A = [A,j];
%    end
% end
% 
% maxD = datasample(A,1);


%update corner node and width/length information for new retrangles
%NC1 = zeros(d);

Width(maxIndex, maxD)=1/3*maxWidth;   
Center = Corner + Width/2;

area(maxIndex)=1;
for j=1:d
    area(maxIndex) = area(maxIndex)*Width(maxIndex, j);
end

if(area(maxIndex)<Vn)
    Vn = area(maxIndex);
end

a = Center(maxIndex,:);
f(maxIndex,1)= cluster(a,d);

if(f(maxIndex,1)<Mn)
    Mn = f(maxIndex,1);
end


if (area(maxIndex) ~= Vn || f(maxIndex) ~= Mn)
    flag = 1;
else 
    flag = 0;
end

NC1 = zeros(1,d);
NC2 = zeros(1,d);
NW1 = zeros(1,d); 
NW2 = zeros(1,d);
NCE1 = zeros(1,d);
NCE2 = zeros(1,d);
BP = zeros(1,d);
NA1= ones(1,1);
NA2= ones(1,1);
NF1= zeros(1,1);
NF2 = zeros(1,1);

for j = 1:d
    if (j == maxD)
        NC1(j) = 1/3*maxWidth + Corner(maxIndex,j);
        NC2(j) = 2/3*maxWidth + Corner(maxIndex,j);
        NW1(j) = 1/3*maxWidth;
        NW2(j) = 1/3*maxWidth;
    else 
        NC1(j) = Corner(maxIndex,j);
        NC2(j) =  Corner(maxIndex,j);
        NW1(j) = Width(maxIndex,j);
        NW2(j) =  Width(maxIndex,j);
    end
    NA1 = NA1*NW1(1,j);
    NA2 = NA2*NW2(1,j);
    %Center(j, maxIndex) = Corner(j,maxIndex)+ 1/2 * Width(j,maxIndex);
    % area(maxIndex+1) = 1/3 * area(maxIndex);
    % area(maxIndex+2) = 1/3 * area(maxIndex);
    % area = area * Width(j,maxIndex);
end
NCE1=NC1+NW1/2;
NCE2=NC2+NW2/2;
NF1=cluster(NCE1,d);
NF2 = cluster(NCE2,d);

if(NF1<Mn)
    Mn = NF1;
    BP = NCE1
    if(NF2<Mn)
        Mn = NF2;
        BP = NCE2
    end
end

if(NA1<Vn)
    Vn = NA1;
    if(NA2<Vn)
        Vn = NA2;
    end
end

%update corner and width information
Corner = [Corner;NC1;NC2]; 
Width = [Width;NW1;NW2];
Center = [Center;NCE1;NCE2];
area = [area;NA1;NA2];
f = [f;NF1;NF2];

[m,d] = size(Corner);

% Vn = min(area);
%  Mn = min(f);


%calculate center nodes and retrangle area
%Center = Corner + Width/2;

%{
for j = 1:m
    for k = 1:d
        Center(j,k)=Corner(j,k)+1/2*Width(j,k);
%        area(j,1)=area(j,1)*Width(j,k);
%        f(j) = funct(Center);
    end
end
%}

%get the smallest function value among all vertices
%{
Mn = intmax('int64');
for j = 1:m
   % for k = 1:d
   %     a(k) = Center(j,k);
   % end
   a = Center(j,:);
        f(j) = funct1(a,d);
    %{   
        if(f(j)<Mn)
            Mn = f(j);
        end
        %}
end
Mn = min(f);
 %}

%M(1,i)=Mn;

% get the volume of smallest rectangle:travsal all rectangles
%{
area = ones(2*n+1,1);
for j = 1:m
    for k = 1:d
        area(j,1)= area(j,1)*Width(j,k);
    end
end


Vn = intmax('int64');
for j = 1 : m
        if (area(j,1) < Vn)
            Vn = area(j,1);
        end
end
Vn;
%}
%V(1,i)=Vn;

% get GVn
%GVn = d*(Vn*log(1/Vn)).^(d/2);
%GVn = d*(Vn*log(1/Vn)).^(d/4);
%GVn = 4*(Vn^(1/2))*(d*log(3/(Vn^(2/d))))^(d/4);

GVn = 4.0*(Vn^(r/d))*(0.5*d*log(evalCount+2.0))^(r/d);

%choose the retrangle which have the largest Rho
%rho = intmin;
%figure (1)
%axis equal 
%axis ([0 1 0 1])
%hold on


if (i==1||flag)
%if (i==1||V(i)~=V(i-1) || M(i)~= M(i-1))
    for j = 1 : m
    % a= Center(j,:);
    % b = Center(j,:);
     %getrho(j) = area(j,1).^(2/d)/(funct(a,d)-M(1,i)+d*GVn);
    %   getrho(j) = area(j,1).^(2/d)/(f(j)-Mn+GVn); %better for 2d raf
  %getrho(j) =area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
   getrho(j) = area(j,1)^(r/d)/(f(j)-Mn+GVn);
   
    end
     
        [tmp, ind] = sort(getrho,'descend');
        Center = Center(ind,:);
        Corner = Corner(ind,:);
        Width = Width(ind,:);
        f = f(ind);
        area = area(ind);
        
      %{  
            if (getrho(j) > rho)
            rho = getrho(j);
            maxIndex = j; 
       %}
        
      %      figure (1);
      %      plot(Center(maxIndex,1),Center(maxIndex,2),'k.');
      %      hold on
      %      Center(maxIndex,:);
            
      %      drawnow; 
       %     pause(0.00000000000001);
        %end

maxIndex = 1;
   
else
    maxIndex = maxIndex+1;
end




 Center(maxIndex, :);
 plot(Center(:,1), Center(:,2),'.');
% scatter3(Center(:,1),Center(:,2),Center(:,3),'*');

% z = f(maxIndex,1);
% x = Center(:,1);
% y = Center(:,2);



end
% scatter3(x,y,z,441,z,'*');
Mn
toc

function y = funct(c,d)
%n = size(c);
sum=0;
for i = 1:d
    sum = sum+ (c(i)-0.5).^2;
end
y = sum;
end

function y = funct1(c,d)
%n = size(c);
nc=zeros(d);
sum=20;
B=4.0;
for i = 1:d
    nc(i)=2*B*c(i)-B-0.1;
    sum = sum+ nc(i).^2-10*cos(2*pi*nc(i));
end
  y=sum;
%y = sum*0.001;
end

function y = cluster(c,d)
global evalCount;
global X;

evalCount = evalCount +1;
sum = 0;

for i = 1:200
    if (X(i,1)-c(1))^2 + (X(i,2)-c(2))^2 < (X(i,1)-c(3))^2 + (X(i,2)-c(4))^2
        sum = sum + (X(i,1)-c(1))^2 + (X(i,2)-c(2))^2;
    else
        sum = sum + (X(i,1)-c(3))^2 + (X(i,2)-c(4))^2;
    end

end

    y = sum;

end


