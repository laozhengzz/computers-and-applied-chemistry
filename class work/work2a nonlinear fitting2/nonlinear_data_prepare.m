
x = 0:0.01:pi+2.8;

y = 3.2.*cos(x) + 2.5.*sin(2*x+1) + 5 + 1.5*rand(1,size(x,2));

results = [x' y'];

scatter(results(:,1),results(:,2))