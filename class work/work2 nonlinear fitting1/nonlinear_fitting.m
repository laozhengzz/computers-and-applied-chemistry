clear all
warning off
run nonlinear_data

beta=[1
    1
    1
    1
];

input = results(:,1);
output = results(:,2);

[b,r,J] = nlinfit(input,output,@logfun,beta);

ra=output-r;


r_square = rsqrgen(ra,output);

RMSD = rmsdgen(ra,output);

polyresult = logfun(b,results(:,1));

scatter(input(:,1),output(:,1))
hold on
plot(input(:,1),polyresult(:,1))


