clear all
x = 0.5:0.5:20;

y = 3.2.*log(x) + 5 + 3*rand(1,size(x,2));

results = [x' y'];