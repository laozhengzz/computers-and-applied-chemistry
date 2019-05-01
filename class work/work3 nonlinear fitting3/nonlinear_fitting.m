clear all
warning off
run nonlinear_data

beta=[1
    1
    1
    1
];

input = results(1:50,1);
output = results(1:50,2);

[b,r,J] = nlinfit(input,output,@polyfun3,beta);

ra=output-r;


r_square = rsqrgen(ra,output);

RMSD = rmsdgen(ra,output);

polyresult = polyfun3(b,results(:,1));


scatter(results(:,1),results(:,2))
hold on
plot(results(:,1),polyresult(:,1))

