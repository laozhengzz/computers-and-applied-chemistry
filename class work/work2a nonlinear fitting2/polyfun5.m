function yhat = polyfun5(beta,x)

b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
b4 = beta(4);
b5 = beta(5);
b6 = beta(6);

x1 = x(:,1);


yhat = (b1.*x1.^5+ b2.*x1.^4+ b3.*x1.^3 + b4.*x1.^2 + b5.*x1 + b6 );

