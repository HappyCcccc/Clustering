d=6;
n=2000;


for j = 1:100
i=0;
xbest = rand(d,1)
fbest = cluster(xbest);
while(i<n)

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 5000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.001;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = rand(d,1); niter = 0; dx = inf;

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = gradient(x);
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;
    newalpha = alpha;
    
    while (any(xnew(:)<0|xnew(:)>1))|(cluster(xnew)>cluster(x))
       newalpha = 0.5*newalpha;
       xnew = x -newalpha*g;
    end
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % plot current point
    %plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end

xopt = x;
fb = cluster(xopt);
if(fb<fbest)
    fbest = fb;
    xbest = xopt;
end
i = i+niter;

end
x = xbest;
toc

cost=cluster(x)
end



function y = gradient(c)

%{
s = 250;
t = 250;
r = 30;
%}
f0 = cluster(c);
g = [0 0 0 0 0 0];
for i =1:6
    cp = c;
    cp(i) = c(i) + 0.01;
    g(i)=(f0-cluster(cp))/0.01;
end

y=g;

end