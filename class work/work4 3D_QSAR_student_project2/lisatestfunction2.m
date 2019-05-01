function yhat = lisatestfunction2(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);



x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);


yhat = (-1/(2.3026*1.9858775*0.001*310))* (b1*x1+ b2*x2+ b3*x3+ b4*x4+ b5);

