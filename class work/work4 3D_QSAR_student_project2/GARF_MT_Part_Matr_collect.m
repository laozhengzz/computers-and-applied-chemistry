%function out = kecsamt_model2(protein_dir,ligand_dir)
clear all
warning off

pKd1=[-7.658107299
-10.90122035
-8.066903062
-7.222058485
-7.521842045
-7.303817638
-3.406631361
-5.069067465
-6.567985264
-5.273465347
-3.583776192
-3.038715174
-5.873032466
-10.51967764
-10.77858163
-10.90122035
-9.702086116
-9.129772047
-9.579447387
-9.920110523
-7.194805434
-10.6559429
-11.90958324
-10.92847341
-8.121409164
-4.469500345
-5.968418144
-9.320543403
-6.540732213
-7.549095096
-7.358323739
-4.387741193
-12.67266866
-13.57201934
-11.71881188
-13.62652544
-10.05637578
-6.704250518
-5.014561363
-9.484061708
-6.404466958
-4.75565738
-8.666470182
-8.175915266
-5.586875432
-9.102518996
-5.723140686
-5.464236703
-6.60886484
-6.813262722
-6.404466958
-5.723140686
-7.208431959
-6.45897306
-6.75875662
-8.175915266
-8.884494589
-13.62652544
-12.05947502
-11.39177527
-12.56365646
-11.95046281
-8.407566198
-4.90554916
-12.61816256
-5.859405941
-13.21772968
-9.129772047
-9.892857472
-8.925374165
-8.843615013
-13.0133318
-9.293290352
-8.31218052
-6.704250518
-8.026023486
-6.731503569
-5.845779415
-5.913912042
-4.387741193
-8.189541791
-12.33200553
-13.21772968
-12.80893392
-11.5689201
-12.1548607
-6.60886484
-11.20100391
-7.658107299
-8.666470182
-9.129772047
-9.334169929
-9.484061708
-6.922274925
-10.84671425
-11.03748561
-12.33200553
-10.43791849
-9.756592217
-4.033451531
-8.829988487
-7.617227723
-9.797471794
-11.60979968
-12.56365646
-11.44628137
-7.140299332
-9.38867603
-8.76185586
-9.293290352
-6.254575178
-10.90122035
-11.10561824
-9.53856781
-9.266037301
-9.53856781
-8.080529588
-10.84671425
-9.865604421
-9.334169929
-8.325807046
-11.09199171
-9.770218743
-8.884494589
-8.993506793
-7.467335943
-7.767119503
-10.90122035
-6.990407552
-5.927538568
-5.723140686
-7.453709417
-6.145562975
-9.53856781
-7.38557679
-10.99660603
-8.829988487
-8.802735436
-8.884494589
-6.1591895
-9.484061708
-3.37937831
-4.960055261
-6.704250518
-11.20100391
-12.40013815
-11.37814875
-8.680096707
-8.925374165
-8.680096707
-7.889758232
-8.175915266
-11.65067925
-12.56365646
-12.56365646
-15.34346765
-12.83618697
-10.84671425
-14.88016578
-15.53423901
-14.93467189
-13.55839282
-12.9588257
-14.03532121
-11.35089569
-13.51751324
-14.74390053
-11.60979968
-8.720976284
-9.020759843
-10.49242459
-8.979880267
-6.118309924
-5.23258577
-7.113046281
-10.51967764
-7.589974672
-9.879230946
-11.85507714
-8.271300944
-6.704250518
-7.767119503
-4.796536956
-8.966253742
-10.28802671
-10.08362883
-9.293290352
-11.41902832
-11.44628137
-7.971517384
-6.963154501
-12.33200553
-9.920110523
-7.903384757
-12.80893392
-8.066903062
-14.15795994
-12.1548607
-6.254575178
-10.69682247
-8.584711029
-10.24714713
-4.333235091
-10.54693069
-10.24714713
-5.450610177
-5.791273313
-10.24714713
-11.85507714
-5.041814414
-9.661206539
-10.28802671
-10.35615934
-10.69682247
-7.440082892
-7.249311536
-8.31218052
-8.802735436
-6.935901451
-10.49242459
-12.80893392
-7.617227723
-11.37814875
-11.90958324
-12.2638729
-15.19357587
-13.21772968
-9.892857472
-14.28059866
-12.38651163
-11.44628137
-11.41902832
-12.2638729
-5.246212296
-5.559622381
-14.33510477
-8.925374165
-7.780746028
-8.216794842
-8.339433571
-6.813262722
-10.6559429
-10.54693069
-12.2093668
-12.11398112
-7.712613401
-7.589974672
-5.178079668
-10.6559429
-8.121409164
-4.891922634
-10.43791849
-6.00929772
-9.266037301
-12.11398112
-11.24188349
-11.05111213
-10.6559429
-6.935901451
-9.702086116
-8.666470182
-10.46517154
-7.453709417
-7.971517384
-7.971517384
-8.325807046
-9.225157725
-7.821625604
-8.121409164
-9.715712641
-9.947363574
-9.225157725
-9.225157725
-9.075265945
-8.748229335
-7.862505181
-7.113046281
-7.821625604
-6.75875662
-10.08362883
-10.49242459
-3.324872208
-7.072166705
-6.949527976
-5.736767212
-5.586875432
-7.521842045
-4.128837209
-7.971517384
-11.85507714
-9.225157725
-6.949527976
-9.334169929
-6.527105687
-7.739866452
-8.448445775
-10.08362883
-9.484061708
-9.865604421
-5.600501957
-9.266037301
-7.848878655
-6.663370942
-5.995671195
-6.268201704
-7.017660603
-4.93280221
-8.543831453
-4.496753396
-7.480962468
-8.393939673
-8.066903062
-10.43791849
-11.65067925
-7.222058485
-9.075265945
-5.750393737
-10.11088188
-7.971517384
-7.440082892
-7.862505181
-9.334169929
-6.486226111
-10.31527976
-6.104683399
-9.225157725
-10.54693069
-14.40323739
-12.90431959
-13.95356205
-7.344697214
-3.202233479
-10.54693069
-8.884494589
-7.821625604
-10.90122035
-8.802735436
-10.08362883
-7.358323739
-10.54693069
-7.072166705
-8.829988487
-7.017660603
-9.38867603
-10.54693069
-12.63178909
-13.72191112
-11.90958324
-7.317444163
-10.16538798
-13.7627907
-7.222058485
-8.175915266
-11.03748561
-5.927538568
-9.83835137
-10.02912273
-9.429555607
-8.734602809
-7.412829841
-10.11088188
-10.46517154
-11.85507714
-11.60979968
-12.44101773
-8.203168317
-9.470435183
-7.017660603
-3.992571955
-4.592139074
-5.450610177
-4.387741193
-5.518742805
-4.619392125
-3.365751784
-11.31001612
-6.75875662
-5.518742805
-7.971517384
-5.750393737
-4.210596362
-8.053276537
-4.15609026
-10.49242459
-6.268201704
-7.930637808
-9.184278149
-8.720976284
-5.23258577
-3.365751784
-4.101584158
-4.75565738
-3.093221276
-10.00186968
-9.293290352
-4.633018651
-6.949527976
-10.28802671
-6.567985264
-4.428620769
-8.434819249
-5.845779415
-8.039650012
-6.445346535
-5.627755008
-4.837416532
-5.68226111
-6.418093484
-4.592139074
-5.464236703
-10.62868985
-7.630854248
-10.28802671
-10.49242459
-10.96935298
-10.90122035
-7.222058485
-10.13813493
-10.46517154
-11.20100391
-11.85507714
-8.61196408
-5.859405941
-7.113046281
-6.45897306
-7.167552383
-11.10561824
-11.99134239
-11.44628137
-11.55529358
-6.949527976
-11.4599079
-8.530204927
-8.584711029
-5.041814414
-8.584711029
-4.087957633
-4.128837209
-9.484061708
-5.396104076
-11.85507714
-7.222058485
-6.567985264
-6.60886484
-8.175915266
-8.421192724
-9.334169929
-9.375049505
-11.60979968
-12.64541561
-9.334169929
-4.196969837
-7.249311536
-7.562721621
-15.94303477
-13.77641722
-8.407566198
-8.925374165
-8.026023486
-7.521842045
-8.066903062
-9.974616624
-11.09199171
-12.56365646
-13.24498273
-5.913912042
-7.276564587
-8.884494589
-9.742965692
-11.95046281
-6.1591895
-5.028187889
-3.025088648
-8.448445775
-12.69992171
-12.05947502
-10.31527976
-12.61816256
-12.65904214
-4.960055261
-8.76185586
-7.821625604
-6.499852636
-10.1245084
-8.598337555
-5.777646788
-3.774547548
-7.658107299
-11.60979968
-8.925374165
-10.45154501
-10.08362883
-9.606700437
-10.96935298
-7.004034078
-6.622491365
-8.026023486
-8.434819249
-8.175915266
-5.041814414
-5.586875432
-11.60979968
-11.65067925
-12.97245222
-9.661206539
-3.815427124
-7.358323739
-6.9086484
-7.262938061
-7.739866452
-11.31001612
-8.734602809
-7.521842045
-11.20100391
-4.837416532
-4.837416532
-5.53236933
-5.859405941
-6.840515773
-9.429555607
-9.048012894
-13.47663366
-9.38867603
-11.25551002
-8.31218052
-9.83835137
-5.600501957
-10.27440018
-7.767119503
-9.38867603
-7.358323739
-7.044913654
-11.14649781
-14.33510477
-11.85507714
-12.50915036
-8.339433571
-8.31218052
-5.396104076
-9.102518996
-9.742965692
-11.85507714
-3.978945429
-3.978945429
-5.041814414
-11.10561824
-7.290191112
-9.429555607
-5.436983652
-12.97245222
-6.567985264
-6.826889247
-6.731503569
-8.271300944
-10.00186968
-10.69682247
-7.480962468
-7.344697214
-6.349960857
-4.128837209
-5.913912042
-5.886658991
-5.178079668
-15.46610638
-10.20626756
-5.98204467
-11.55529358
-11.32364264
-7.113046281
-15.53423901
-3.324872208
-15.91578172
-8.393939673
-15.76588994
-15.2617085
-15.79314299
-15.97028782
-7.358323739
-12.11398112
-12.2638729
-11.5689201
-11.5689201
-7.426456367
-9.674833065
-7.562721621
-9.102518996
-8.175915266
-9.98824315
-10.6559429
-7.371950265
-8.993506793
-12.2638729
-6.663370942
-9.974616624
-9.53856781
-9.797471794
-6.404466958
-11.17375086
-7.140299332
-11.85507714
-10.58781027
-6.254575178
-11.4054018
-11.4054018
-11.4054018
-11.4054018
-6.527105687
-6.377213907
-8.884494589
-6.540732213
-7.480962468
-8.720976284
-13.7627907
-14.25334561
-12.80893392
-14.32147824
-11.54166705
-6.104683399
-7.68536035
-11.85507714
-8.530204927
-9.293290352
-9.129772047
-8.884494589
-7.222058485
-6.813262722
-13.7627907
-12.2638729
-8.31218052
-12.76805434
-11.85507714
-7.821625604
-14.58038222
-14.83928621
-9.920110523
-7.726239926
-10.38341239
-9.865604421
-6.036550771
-8.979880267
-11.03748561
-10.38341239
-8.952627216
-5.995671195
-6.567985264
-8.979880267
-10.31527976
-12.27749942
-11.84145061
-15.39797375
-9.783845268
-7.903384757
-8.584711029
-10.13813493
-8.625590606
-9.102518996
-8.434819249
-11.36452222
-11.82782408
-12.40013815
-12.35925858
-10.84671425
-6.404466958
-7.794372554
-7.358323739
-8.421192724
-10.19264103
-9.920110523
-10.05637578
-9.947363574
-13.89905595
-14.22609256
-14.78478011
-14.40323739
-13.91268248
-10.58781027
-12.46827078
-9.742965692
-8.325807046
-5.804899839
-5.450610177
-6.854142298
-6.922274925
-11.60979968
-9.98824315
-6.091056873
-10.13813493
-7.971517384
-9.38867603
-8.175915266
-3.542896615
-11.55529358
-12.75442781
-12.93157265
-10.62868985
-9.252410776
-4.387741193
-12.3728851
-11.02385908
-4.660271702
-13.42212756
-5.873032466
-11.85507714
-11.55529358
-9.429555607
-4.987308312
-4.087957633
-11.31001612
-5.055440939
-13.21772968
-11.90958324
-10.90122035
-11.60979968
-3.733667971
-6.976781027
-5.396104076
-8.31218052
-12.75442781
-9.020759843
-11.20100391
-7.358323739
-10.28802671
-11.95046281
-12.63178909
-11.71881188
-12.05947502
-11.44628137
-7.044913654
-10.38341239
-7.589974672
-8.543831453
-8.244047893
-6.840515773
-10.84671425
-8.530204927
-3.67916187
-11.85507714
-11.10561824
-4.224222887
-4.987308312
-7.480962468
-6.75875662
-5.491489754
-4.564886023
-4.237849413
-3.951692379
-9.811098319
-10.30165324
-10.72407552
-9.933737048
-10.31527976
-10.7649551
-10.68319595
-10.32890629
-9.565820861
-6.676997467
-6.60886484
-8.870868064
-8.639217131
-8.652843656
-5.028187889
-14.55312917
-10.58781027
-10.90122035
-8.91174764
-11.80057103
-11.20100391
-11.99134239
-11.31001612
-10.79220815
-11.39177527
-9.266037301
-14.43049044
-11.21463044
-7.971517384
-9.68845959
-10.90122035
-12.05947502
-8.584711029
-9.129772047
-10.24714713
-11.85507714
-10.46517154
-9.034386369
-7.930637808
-9.157025098
-4.90554916
-6.404466958
-5.886658991
-8.203168317
-14.29422519
-4.442247294
-3.22948653
-11.77331798
-14.60763528
];

%%% In the entire code, number 1 associated with anything reflects the
%%% sampling in the region and number 2 means the sampling at the local
%%% minima of the region in number 1.

%function GARF_MT_EXE(home_dir)
home_dir=('C:\Users\laozh\Documents\test_set_795_PL\');

%vdw_list_num_AF = contact_list_final_AF_format(vdw_all_AF_contact_list);

run GARF_Potential
%run GARF_Potential_mute_nitrogen
%run GARF_nonbonding_PP

%run Solvation_Mat
%run GARF_inter_PP

%run LL_torsion_FF

rd = 0.005;

R_srch=0.5;
R_L_tor_srch=0.1;
hnstate = 0.5/rd;
halva=hnstate;
halv=100;
halv1=10;
RcutoffPL = 6;
const=halva*2;
const1=halv1*2;
column_num=100;

size_vdw_d=size(vdw_d);
for i=2:size_vdw_d(1,2)
    vdw_d(vdw_d(:,i)<=10^-10,i)=10^-10;
end

C=1  ;
Cs=2 ;
CA=3 ;
CB=4 ;
CC=5 ;
CN=6 ;
CR=7 ;
CT=8 ;
CV=9 ;
CW=10;
H=11 ;
HO=12;
N=13;
N2=14;
N3=15;
NA=16;
NB=17;
O=18 ;
O2=19;
OH=20;
S=21;
SH=22;



C_3=1;
C_2=2;
C_1=3;
C_ar=4;
O_3=5;
O_3p=6;
O_2=7;
O_co2=8;
O_2v=9;
N_2=10;
N_am=11;
N_pl3=12;
N_4=13;
P=14;
F=15;
Cl=16;
Br=17;
I=18;
C_cat=19;
S_3=20;
S_o=21;
HNa=22;
HOa=23;
HSa=24;
H3 =25;
C_3Oa=26;
C_3Na=27;
C_3La=28;
C_2Oa=29;
C_2Na=30;
C_2La=31;
C_arOa=32;
C_arNa=33;
C_arLa=34;
HCa=35;

N_3=12;
N_1=10;
N_ar=10;
S_2=21;
S_o2=22;





K    = 29 ;
MG   = 30  ;
CAA   = 31	;
AL   = 32	;
MN   = 33  	;
FE   = 34    ;
CO   = 35    ;
NI   = 36    ;
CU   = 37    ;
ZN   = 38    ;


N_A=0;
D=1;
D2=2;
A=3;
DA=4;

ALA=101;
ARG=102;
ASN=103;
ASP=104;
CYS=105;
GLN=106;
GLU=107;
GLY=108;
HIS=109;
ILE=110;
LEU=111;
LYS=112;
MET=113;
PHE=114;
PRO=115;
SER=116;
THR=117;
TRP=118;
TYR=119;
VAL=120;



carbon='C';
oxygen='O';
nitrogen='N';
sulfer='S';
fluorin='F';
chlorine='Cl';
bromine='Br';
iodine='I';
phosphor='P';
hydrogen='H';

count=1

%APO_protein_dir1 = 'C:\Users\John Zheng\Dropbox\Test_cases\Streptavidin_apo_holo\Apo_protein\';
%list_ApoP=dir(HOLO_protein_dir1);

%home_dir = 'F:\Streptavidin_apo_holo\Test_set\';
list=dir(home_dir);
size_list=size(list);
%Aout=zeros(200,6);
%Bout=zeros(200,28);

tn=0;
clear test_name
%for zh=15:15

%fp1=fopen(strcat(home_dir,'PL_MT_binding.txt'),'wt');
%fp2=fopen(strcat(home_dir,'P_MT.txt'),'wt');
%fp3=fopen(strcat(home_dir,'PL_MT.txt'),'wt');

for zh=3:size_list(1,1)
    
    tic
    if length(list(zh,1).name)<8
        continue
    elseif length(list(zh,1).name)>=8
        list_name=list(zh,1).name;
        
        if strcmp(list_name(end-3:end),'.pdb')==1
            protein_dir=[ strcat(home_dir,list_name)];
            
            id_name = list_name(1:4);
            ligand_pot = strcat(id_name,'_ligand.mol2');
            ligand_dir=[ strcat(home_dir,ligand_pot)];
            
            [ligand,name_str] = ligand_define_sol(ligand_dir);
            
            [protein,protein_Namelist,protein_sol] = protein_define_final_sol(protein_dir);
            
            protein(protein(:,5)==121,:)=[];
            
            [sum_com_dp_dp,sum_com_dp_indp,sum_com_indp_indp,sum_com_sc_sc] = contactdetectGARF(protein,ligand,RcutoffPL,vdw_d);
            
            
            
            sum_com_dp_dp_all(count,1) = sum_com_dp_dp;
            sum_com_dp_indp_all(count,1) = sum_com_dp_indp;
            sum_com_indp_indp_all(count,1) = sum_com_indp_indp;
            sum_com_sc_sc_all(count,1) = sum_com_sc_sc;
            
            
            
            toc
            %fprintf(fp,' % s  % 5d  % 5d\n',ligand_pot,Aout(count,1), Bout(count,1));
            count=count+1;
            
            
        end
    end
end

beta=[0.1
    0.1
0.1
0.1
0.1

];




train1 = [sum_com_dp_dp_all sum_com_dp_indp_all sum_com_indp_indp_all sum_com_sc_sc_all];


options.MaxIter = 10000;
options.TolFun = 1e-10;
options.TolX = 1e-10;
[b,r,J] = nlinfit(train1,pKd1,@lisatestfunction2,beta,options);

c=nlparci(b,r,J);
ra=pKd1-r;


size_r=size(r);

for i=1:1:size_r(1,1)
    tmp_simu(i)=ra(i)-mean(ra);
    tmp_pKd(i)=pKd1(i)-mean(pKd1);
    fenzi(i)=tmp_simu(i)*tmp_pKd(i);
    fenmu1(i)=tmp_simu(i)*tmp_simu(i);
    fenmu2(i)=tmp_pKd(i)*tmp_pKd(i);
end


rs=sum(fenzi)/(sqrt(sum(fenmu1))*sqrt(sum(fenmu2)));

r_square=rs^2;


for i=1:1:size_r(1,1)
    ss(i)=(ra(i)-pKd1(i))^2;
end

RMSD=sqrt(sum(ss)/795);

















