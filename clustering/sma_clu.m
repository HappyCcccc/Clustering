

for i = 1:100
    
ObjectiveFunction = @cluster;
%X0 = [0.5 0.5 0.5 0.5 0.5 0.5];% Starting point
X0 = rand(6,1)
lb =[0 0 0 0 0 0];
ub = [1 1 1 1 1 1];
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0,lb,ub)

end