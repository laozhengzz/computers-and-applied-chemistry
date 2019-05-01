function yhat = polyfun7(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);
b6 = beta(6);
b7 = beta(7);
b8 = beta(8);

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
x5 = x(:,5);
x6 = x(:,6);
x7 = x(:,7);

yhat = (b1*x1^7+ b2*x2^6+ b3*x3^5 + b4*x4^4 + b5*x5^3 + b6*x6^2 + b7*x7 + b8 );

