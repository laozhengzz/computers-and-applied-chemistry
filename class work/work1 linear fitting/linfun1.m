function yhat = linfun1(beta,x)

b1 = beta(1);
b2 = beta(2);



x1 = x(:,1);


yhat = (b1*x1+ b2);

