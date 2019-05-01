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

yhat = (b1.*x1.^7+ b2.*x1.^6+ b3.*x1.^5 + b4.*x1.^4 + b5.*x1.^3 + b6.*x1.^2 + b7.*x1 + b8 );

