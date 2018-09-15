x1 = zeros(10000);
y1 = zeros(10000);
c3 = zeros(1,1);
k=0;
figure
for i=1:100
    for j=1:100
       
        
        x1(k)= i*0.01
        y1(k)= j*0.01
  %      c2(i*j,1) = j*0.01
       c3(i*j,1) = cluster([0.2,0.2,i*0.01,j*0.01]);
    k=k+1;
    end
end
%plot3(c1,c2,c3,'*');

function y = cluster(c)
%global evalCount;
global X;

%evalCount = evalCount +1;
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
