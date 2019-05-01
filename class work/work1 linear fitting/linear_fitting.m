clear all
warning off
run linear_data

beta=[1
    1
];

input = results(:,1);
output = results(:,2);

[b,r,J] = nlinfit(input,output,@linfun1,beta);

ra=output-r;

r_square = rsqrgen(ra,output);

RMSD = rmsdgen(ra,output);