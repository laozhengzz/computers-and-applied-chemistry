function yhat = polyfun3(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);


yhat = (b1*x1^3+ b2*x2^2+ b3*x3 + b4 );

