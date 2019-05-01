function yhat = polyfun3(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);

x1 = x(:,1);

yhat = (b1.*x1.^3+ b2.*x1.^2+ b3.*x1 + b4 );

