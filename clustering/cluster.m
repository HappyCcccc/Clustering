function y = cluster(c)
%global evalCount;
global Z;
global k;
global d1;
global datan;

%evalCount = evalCount +1;
sum = 0;
datan=100;
k=5;
d1=2;
for i = 1:datan
    mindist = intmax('int64');
    for j = 1:k
        dist = 0;
        
        for l=1:d1
            dist = dist+(c((j-1)*d1+l)-Z(i,l))^2; 
       %     sqrt
            
        end
        
        if dist < mindist
            mindist = dist;
        end
    end
   sum = sum + mindist; 
end
%sum = sqrt(sum);
y = sum;

end
