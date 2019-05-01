function yhat = polyfun5(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);
b6 = beta(6);

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
x5 = x(:,5);


yhat = (b1*x1^5+ b2*x2^4+ b3*x3^3 + b4*x4^2 + b5*x5 + b6 );

