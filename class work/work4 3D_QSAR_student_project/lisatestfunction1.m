function yhat = lisatestfunction1(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);
b6 = beta(6);
b7 = beta(7);
b8 = beta(8);
b9 = beta(9);



x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
x5 = x(:,5);
x6 = x(:,6);
x7 = x(:,7);
x8 = x(:,8);


yhat = (-1/(2.3026*1.9858775*0.001*310))* (b1*x1+ b2*x2+ b3*x3+ b4*x4+ b5*x5+ b6*x6+ b7*x7+ b8*x8 + b9);

